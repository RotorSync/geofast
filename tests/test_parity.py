"""Parity check: geofast's Rust-backed generator vs the upstream field-predictor
`quote` binary on the same field.

Both call the identical `ingest::to_local_feet` + `quote()` with default params, so
the chosen angle and total spray distance must match to within rounding. The test
is skipped unless FP_QUOTE_BIN points at a built upstream `quote` binary.

Run:
    FP_QUOTE_BIN=~/Desktop/field-predictor/rust/target/release/quote pytest tests/test_parity.py
"""
import os
import re
import subprocess

import pytest

from geofast.spray_line_generator import SprayLineGenerator, SprayConfig

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")
QUOTE_BIN = os.environ.get("FP_QUOTE_BIN")

# Same field as fixtures/field_irregular.kml
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


@pytest.mark.skipif(not QUOTE_BIN or not os.path.exists(QUOTE_BIN),
                    reason="FP_QUOTE_BIN not set / binary not built")
def test_parity_with_upstream_quote():
    kml = os.path.join(FIXTURES, "field_irregular.kml")
    out = subprocess.run([QUOTE_BIN, kml], capture_output=True, text=True, check=True).stdout

    # Parse per-block "spray=NNNN ft" and "angle=NN.N°"
    sprays = [float(x) for x in re.findall(r"spray=(\d+)\s*ft", out)]
    angles = [float(x) for x in re.findall(r"angle=([-\d.]+)°", out)]
    assert sprays, f"could not parse quote output:\n{out}"
    upstream_spray_ft = sum(sprays)

    res = SprayLineGenerator(SprayConfig(swath_width_ft=50.0)).generate(L_FIELD)

    # spray_ft is summed across blocks; quote prints %.0f so allow 1 ft/block rounding
    assert res.total_spray_distance_ft == pytest.approx(upstream_spray_ft, abs=len(sprays) + 1)
    # dominant angle matches the first block
    assert res.spray_bearing_deg == pytest.approx(angles[0], abs=0.05)
