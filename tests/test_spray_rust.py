"""Tests for the Rust-backed spray line generator (geofast._native).

These exercise the seam (SprayLineGenerator.generate), the optimizer path, and
the end-to-end formats API. Parity against the upstream field-predictor `quote`
binary lives in test_parity.py (skipped if the binary isn't available).
"""
import json
import math
import os

import pytest

from geofast.spray_line_generator import SprayLineGenerator, SprayConfig

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")

# L-shaped field (lon, lat), matches fixtures/field_irregular.geojson
L_FIELD = [
    [
        (-97.500, 38.000),
        (-97.488, 38.000),
        (-97.488, 38.004),
        (-97.494, 38.004),
        (-97.494, 38.009),
        (-97.500, 38.009),
        (-97.500, 38.000),
    ]
]

# A simple rectangle ~2878 ft (E-W) x ~2550 ft (N-S) => ~168 acres
RECT_FIELD = [
    [
        (-97.500, 38.000),
        (-97.490, 38.000),
        (-97.490, 38.007),
        (-97.500, 38.007),
        (-97.500, 38.000),
    ]
]


def test_native_module_imports():
    import geofast._native as nat
    assert hasattr(nat, "plan_lines")


def test_generate_basic_rectangle():
    gen = SprayLineGenerator(SprayConfig(swath_width_ft=50.0))
    res = gen.generate(RECT_FIELD)

    assert res.num_lines > 0
    assert len(res.lines) == res.num_lines
    assert res.total_spray_distance_ft > 0
    assert res.total_spray_distance_miles == pytest.approx(res.total_spray_distance_ft / 5280.0)
    # Analytic area ~168 acres
    assert res.field_area_acres == pytest.approx(168.5, rel=0.02)
    # Each spray line is a 2-point (lon,lat) segment
    for ln in res.lines:
        assert len(ln) == 2
        assert len(ln[0]) == 2 and len(ln[1]) == 2


def test_generate_irregular_field():
    gen = SprayLineGenerator(SprayConfig(swath_width_ft=50.0))
    res = gen.generate(L_FIELD)
    assert res.num_lines > 0
    assert res.total_spray_distance_ft > 0
    assert res.field_area_acres > 0


def test_bearing_override_is_ignored():
    """The engine owns angle selection; bearing_override must not change output."""
    gen = SprayLineGenerator(SprayConfig(swath_width_ft=50.0))
    a = gen.generate(RECT_FIELD, bearing_override=0)
    b = gen.generate(RECT_FIELD, bearing_override=90)
    assert a.num_lines == b.num_lines
    assert a.total_spray_distance_ft == pytest.approx(b.total_spray_distance_ft)
    assert a.spray_bearing_deg == pytest.approx(b.spray_bearing_deg)


def test_determinism():
    gen = SprayLineGenerator(SprayConfig(swath_width_ft=50.0))
    r1 = gen.generate(L_FIELD)
    r2 = gen.generate(L_FIELD)
    assert r1.lines == r2.lines
    assert r1.total_spray_distance_ft == r2.total_spray_distance_ft


def test_swath_affects_line_count():
    """Wider swath => fewer passes."""
    narrow = SprayLineGenerator(SprayConfig(swath_width_ft=30.0)).generate(RECT_FIELD)
    wide = SprayLineGenerator(SprayConfig(swath_width_ft=80.0)).generate(RECT_FIELD)
    assert wide.num_lines < narrow.num_lines


def test_lines_lie_within_bbox():
    """Generated spray lines stay within the field bounding box (small epsilon)."""
    gen = SprayLineGenerator(SprayConfig(swath_width_ft=50.0))
    res = gen.generate(RECT_FIELD)
    lons = [c[0] for ring in RECT_FIELD for c in ring]
    lats = [c[1] for ring in RECT_FIELD for c in ring]
    minlon, maxlon, minlat, maxlat = min(lons), max(lons), min(lats), max(lats)
    eps = 1e-4
    for ln in res.lines:
        for (lon, lat) in ln:
            assert minlon - eps <= lon <= maxlon + eps
            assert minlat - eps <= lat <= maxlat + eps


# Rectangle with a small interior hole (~460x365 ft, well under the ~1760 ft
# turn-around break-even) so scan lines fly straight across it -> transit "hops".
HOLED_FIELD = [
    [
        (-97.500, 38.000),
        (-97.490, 38.000),
        (-97.490, 38.008),
        (-97.500, 38.008),
        (-97.500, 38.000),
    ],
    [
        (-97.4958, 38.0035),
        (-97.4942, 38.0035),
        (-97.4942, 38.0045),
        (-97.4958, 38.0045),
        (-97.4958, 38.0035),
    ],
]


def test_transit_hops_exposed():
    """A field with a small internal gap should expose boom-off transit hops."""
    res = SprayLineGenerator(SprayConfig(swath_width_ft=50.0)).generate(HOLED_FIELD)
    assert res.num_lines > 0
    assert len(res.transit_lines) > 0, "expected fly-through hops over the hole"
    assert res.transit_distance_ft > 0
    for ln in res.transit_lines:
        assert len(ln) == 2 and len(ln[0]) == 2


def test_solid_field_has_no_hops():
    """A simple convex field needs no fly-through hops."""
    res = SprayLineGenerator(SprayConfig(swath_width_ft=50.0)).generate(RECT_FIELD)
    assert res.transit_lines == []
    assert res.transit_distance_ft == 0.0


def test_geojson_emits_hop_features():
    from geofast.spray_optimizer import generate_spray_pattern_geojson
    geom = {"type": "Polygon", "coordinates": HOLED_FIELD}
    fc = generate_spray_pattern_geojson(geom, config={"swath_width_ft": 50.0})
    hops = [f for f in fc["features"] if f["properties"].get("type") == "hop"]
    assert len(hops) > 0
    assert fc["properties"]["num_hops"] == len(hops)
    assert fc["properties"]["hop_feet"] > 0
    for f in hops:
        assert f["geometry"]["type"] == "LineString"


def test_end_to_end_generate_spray_patterns(tmp_path):
    from geofast.formats import generate_spray_patterns

    out = tmp_path / "out.geojson"
    generate_spray_patterns(
        os.path.join(FIXTURES, "field_irregular.geojson"),
        str(out),
        config={"swath_width_ft": 50.0},
    )
    data = json.loads(out.read_text())
    assert data["type"] == "FeatureCollection"
    spray = [f for f in data["features"]
             if f.get("properties", {}).get("feature_type") == "spray_line"
             or f.get("properties", {}).get("type") == "spray_line"]
    assert len(spray) > 0
    for f in spray:
        assert f["geometry"]["type"] == "LineString"
