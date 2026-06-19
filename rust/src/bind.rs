//! PyO3 binding: expose the field-predictor planner to geofast's Python layer.
//!
//! The geometry/planning algorithm is the unchanged field-predictor crate; this
//! module is only the lon/lat <-> feet glue and the Python entry point. The seam
//! it serves is `geofast.spray_line_generator.SprayLineGenerator.generate`.

use pyo3::prelude::*;

use crate::geom::{Geom, Polygon, Pt};
use crate::ingest::{local_feet_to_lonlat, lonlat_to_local_feet, to_local_feet};
use crate::params::{AircraftParams, ACRE_FT2};
use crate::quote;

/// GeoJSON rings repeat the first vertex as the last; the crate stores rings
/// unclosed (the `coords[:-1]` convention, see ingest::close_to_open). Drop a
/// trailing duplicate so we match what the planner expects.
fn ring_from_coords(coords: Vec<(f64, f64)>) -> Vec<Pt> {
    let mut r: Vec<Pt> = coords.into_iter().map(|(lon, lat)| [lon, lat]).collect();
    if r.len() >= 2 && r[0] == r[r.len() - 1] {
        r.pop();
    }
    r
}

/// Plan spray lines for one field.
///
/// Inputs are WGS84 lon/lat degrees (GeoJSON order: lon first). The crate's own
/// ellipsoidal aeqd projection (`to_local_feet`) handles lon/lat -> feet, runs
/// the full planner (angle search + boustrophedon decomposition + simulate), and
/// the result lines are un-projected back to lon/lat.
///
/// `swath_ft` <= 0 falls back to the param default (50 ft). `trailer_lonlat` and
/// `gpa` only affect block ordering / load timing, not the geometric pass lines,
/// so they default sensibly when omitted.
///
/// Returns `(lines, total_spray_distance_ft, field_area_acres, bearing_deg,
/// total_transit_distance_ft, num_turns, num_runs)`, where each line is
/// `(lon0, lat0, lon1, lat1, is_spray)`. `is_spray == false` segments are the
/// boom-off transit "hops" the planner chose to fly straight across (gaps
/// shorter than the turn-around break-even); `total_transit_distance_ft` is
/// their summed length in feet. `num_turns` is the engine's simulated
/// turn-around count and `num_runs` the total spray passes across plan blocks.
#[pyfunction]
#[pyo3(signature = (exterior, holes, swath_ft, trailer_lonlat=None, gpa=None))]
fn plan_lines(
    exterior: Vec<(f64, f64)>,
    holes: Vec<Vec<(f64, f64)>>,
    swath_ft: f64,
    trailer_lonlat: Option<(f64, f64)>,
    gpa: Option<f64>,
) -> PyResult<(Vec<(f64, f64, f64, f64, bool)>, f64, f64, f64, f64, i64, i64)> {
    let polygon = Polygon {
        exterior: ring_from_coords(exterior),
        interiors: holes.into_iter().map(ring_from_coords).collect(),
    };
    let geom_ll = Geom::Poly(polygon);
    let (geom_ft, center) = to_local_feet(&geom_ll);

    let mut p = AircraftParams::default();
    if swath_ft > 0.0 {
        p.swath_ft = swath_ft;
    }
    if let Some(g) = gpa {
        if g > 0.0 {
            p.gpa = g;
        }
    }

    // Trailer in the same projected frame; default to the projection origin.
    let trailer_ft: Pt = match trailer_lonlat {
        Some((lon, lat)) => lonlat_to_local_feet([lon, lat], center),
        None => [0.0, 0.0],
    };

    let (plans, sim, _t) = quote(&geom_ft, &p, Some(trailer_ft), 2.0);

    let total_spray_distance_ft: f64 = plans.iter().map(|pl| pl.spray_ft).sum();
    let total_transit_distance_ft: f64 = plans.iter().map(|pl| pl.dead_ft).sum();
    let field_area_acres = geom_ft.area() / ACRE_FT2;
    let bearing_deg = plans.first().map(|pl| pl.angle_deg).unwrap_or(0.0);
    // Engine's own counts: total turn-arounds from the flight sim, and spray runs
    // (passes) across all plan blocks. count_effective_lines in the Python layer
    // mis-counts the engine's segmented output, so callers should prefer these.
    let num_turns: i64 = sim.n_turns;
    let num_runs: i64 = plans.iter().map(|pl| pl.n_runs).sum();

    let mut lines: Vec<(f64, f64, f64, f64, bool)> = Vec::new();
    for pl in &plans {
        for line in &pl.lines {
            for &(a, b, is_spray, _l) in line {
                let la = local_feet_to_lonlat(a, center);
                let lb = local_feet_to_lonlat(b, center);
                lines.push((la[0], la[1], lb[0], lb[1], is_spray));
            }
        }
    }

    Ok((
        lines,
        total_spray_distance_ft,
        field_area_acres,
        bearing_deg,
        total_transit_distance_ft,
        num_turns,
        num_runs,
    ))
}

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(plan_lines, m)?)?;
    Ok(())
}
