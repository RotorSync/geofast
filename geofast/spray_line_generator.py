"""
Spray Line Generator for Aerial Crop Dusting

Generates optimal spray lines for agricultural fields using the field-predictor
planning engine (a Rust extension, ``geofast._native``). The engine performs a
proper ellipsoidal projection, a sweep-angle search, per-gap fly-through-vs-turn
decisions, and boustrophedon cell decomposition — the same algorithm used in
RotorSync production.

Coordinates are WGS84 lon/lat degrees throughout (GeoJSON order: lon first).
"""

import json
from typing import List, Tuple, Optional
from dataclasses import dataclass

from . import _native


@dataclass
class SprayConfig:
    """Configuration for spray line generation.

    Only ``swath_width_ft`` is consumed by the planning engine; the remaining
    fields are retained for backward compatibility with existing callers (the
    engine owns angle selection, edge buffering, and short-line handling).
    """
    swath_width_ft: float = 50.0  # Width of spray coverage per pass
    buffer_ft: float = 0.0  # Buffer inside field boundary (0 = spray to edge)
    prefer_cardinal: bool = True  # Prefer N-S or E-W when close
    cardinal_tolerance_deg: float = 15.0  # Snap to cardinal if within this
    min_line_length_ft: float = 50.0  # Skip lines shorter than this


@dataclass
class SprayResult:
    """Result of spray line generation"""
    lines: List[List[Tuple[float, float]]]  # List of spray lines (lon, lat coords)
    total_spray_distance_ft: float
    total_spray_distance_miles: float
    num_lines: int
    field_area_acres: float
    efficiency_miles_per_acre: float
    spray_bearing_deg: float
    estimated_swath_width_ft: float


class SprayLineGenerator:
    """Generates optimal spray lines for agricultural fields.

    Thin wrapper over the ``geofast._native`` planning engine; the geometry and
    planning are performed in Rust.
    """

    def __init__(self, config: Optional[SprayConfig] = None):
        self.config = config or SprayConfig()

    def generate(self, geojson_coords: List, bearing_override: Optional[float] = None) -> SprayResult:
        """
        Generate spray lines for a field.

        Args:
            geojson_coords: GeoJSON polygon coordinates [[outer_ring], [hole1], ...]
                            in WGS84 lon/lat degrees.
            bearing_override: Accepted for API compatibility but ignored — the
                              planning engine selects the spray angle (and per-cell
                              angles) itself via its sweep-angle search.

        Returns:
            SprayResult with spray lines and metrics.
        """
        rings = list(geojson_coords)
        exterior = [(c[0], c[1]) for c in rings[0]]
        holes = [[(c[0], c[1]) for c in ring] for ring in rings[1:]]

        segments, total_distance_ft, area_acres, bearing = _native.plan_lines(
            exterior, holes, self.config.swath_width_ft
        )

        # Each spray segment (is_spray=True) becomes a 2-point line; boom-off
        # transit segments are dropped from the spray-line set.
        line_coords: List[List[Tuple[float, float]]] = [
            [(lon0, lat0), (lon1, lat1)]
            for (lon0, lat0, lon1, lat1, is_spray) in segments
            if is_spray
        ]

        total_distance_miles = total_distance_ft / 5280.0
        efficiency = total_distance_miles / area_acres if area_acres > 0 else 0.0

        return SprayResult(
            lines=line_coords,
            total_spray_distance_ft=total_distance_ft,
            total_spray_distance_miles=total_distance_miles,
            num_lines=len(line_coords),
            field_area_acres=area_acres,
            efficiency_miles_per_acre=efficiency,
            spray_bearing_deg=bearing,
            estimated_swath_width_ft=self.config.swath_width_ft
        )

    def generate_geojson(self, geojson_coords: List, bearing_override: Optional[float] = None) -> dict:
        """
        Generate spray lines and return as GeoJSON FeatureCollection.
        """
        result = self.generate(geojson_coords, bearing_override)

        features = []

        # Add original field boundary
        features.append({
            "type": "Feature",
            "properties": {
                "type": "field_boundary",
                "area_acres": result.field_area_acres
            },
            "geometry": {
                "type": "Polygon",
                "coordinates": geojson_coords
            }
        })

        # Add spray lines
        for i, line_coords in enumerate(result.lines):
            features.append({
                "type": "Feature",
                "properties": {
                    "type": "spray_line",
                    "line_number": i + 1,
                    "bearing_deg": result.spray_bearing_deg
                },
                "geometry": {
                    "type": "LineString",
                    "coordinates": line_coords
                }
            })

        return {
            "type": "FeatureCollection",
            "properties": {
                "total_spray_distance_miles": result.total_spray_distance_miles,
                "num_lines": result.num_lines,
                "field_area_acres": result.field_area_acres,
                "efficiency_miles_per_acre": result.efficiency_miles_per_acre,
                "spray_bearing_deg": result.spray_bearing_deg,
                "swath_width_ft": result.estimated_swath_width_ft
            },
            "features": features
        }


def process_field_file(input_path: str, output_path: str, swath_width: float = 50.0):
    """
    Process a GeoJSON file with field boundaries and generate spray lines.
    """
    config = SprayConfig(swath_width_ft=swath_width)
    generator = SprayLineGenerator(config)

    with open(input_path) as f:
        data = json.load(f)

    results = []

    for feat in data['features']:
        if feat['geometry']['type'] != 'Polygon':
            continue

        coords = feat['geometry']['coordinates']
        props = feat['properties']

        try:
            result = generator.generate_geojson(coords)

            # Merge original properties
            result['properties']['jobId'] = props.get('jobId')
            result['properties']['name'] = props.get('name')
            result['properties']['address'] = props.get('address')
            result['properties']['original_area'] = props.get('area')

            results.append(result)

        except Exception as e:
            print(f"Error processing job {props.get('jobId')}: {e}")

    # Combine all results
    all_features = []
    for r in results:
        all_features.extend(r['features'])

    output = {
        "type": "FeatureCollection",
        "features": all_features
    }

    with open(output_path, 'w') as f:
        json.dump(output, f)

    print(f"Generated spray lines for {len(results)} fields")
    print(f"Output saved to {output_path}")


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("Usage: python spray_line_generator.py <input.geojson> <output.geojson> [swath_width_ft]")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    swath_width = float(sys.argv[3]) if len(sys.argv) > 3 else 50.0

    process_field_file(input_path, output_path, swath_width)
