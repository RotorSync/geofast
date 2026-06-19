//! Top-level planning: cost aggregation + plan_field (joint vs per-part BCD).

use crate::decompose::bcd_plan;
use crate::geom::{poly_distance, Geom, Polygon};
use crate::params::{AircraftParams, ACRE_FT2};
use crate::plan::PlanResult;
use crate::scanline::best_single_direction;

/// `_block_cost`: air + turn time for one block (loads excluded — invariant
/// across decompositions). Each turn adds clock time (t_turn_s) AND the
/// outside-field headland arc traversal (arc_ft_per_turn / ferry_fps). Breaks
/// cost 2 turns; n_runs already counts 1.
pub fn block_cost(plan: &PlanResult, p: &AircraftParams) -> f64 {
    let n_turns = (plan.n_runs - 1).max(0) + plan.n_breaks;
    plan.spray_ft / p.speed_fps
        + plan.dead_ft / p.dead_fps
        + n_turns as f64 * (p.t_turn_s + p.arc_ft_per_turn / p.ferry_fps)
}

/// `combined_time`: total job time for a set of cell plans sharing loads/ferry.
pub fn combined_time(plans: &[PlanResult], p: &AircraftParams, extra_transit_ft: f64) -> f64 {
    let area: f64 = plans.iter().map(|r| r.area_ac).sum();
    let n_loads = if area > 0.0 {
        (area * p.gpa / p.tank_gal).ceil()
    } else {
        0.0
    };
    plans.iter().map(|r| block_cost(r, p)).sum::<f64>()
        + n_loads * (p.t_load_s + p.t_cycle_s)
        + 2.0 * n_loads * p.ferry_ft / p.ferry_fps
        + (plans.len() as f64 - 1.0) * p.t_turn_s
        + extra_transit_ft / p.ferry_fps
}

/// `plan_field`: simplify, then compare joint single-direction vs per-part
/// BCD+merge with inter-part ferry. Returns (plans, total_time_s).
pub fn plan_field(geom: &Geom, p: &AircraftParams, step_deg: f64) -> (Vec<PlanResult>, f64) {
    // TODO: topology-preserving Douglas-Peucker simplify(10.0). Plain DP is
    // equivalent at this tolerance on well-formed fields; validate per
    // PORTING.md §5. The scaffold skips simplification (identity) to avoid a
    // silent vertex divergence — wire in DP once parity is being checked.
    let simplified = geom.clone();

    let joint = best_single_direction(&simplified, p, step_deg, 8);
    let mut best_plans = vec![joint.clone()];
    let mut best_t = combined_time(&best_plans, p, 0.0);

    let mut parts: Vec<Polygon> = simplified
        .parts()
        .into_iter()
        .filter(|g| g.area() > 0.25 * ACRE_FT2)
        .cloned()
        .collect();
    parts.sort_by(|a, b| b.area().partial_cmp(&a.area()).unwrap());

    let mut all_plans: Vec<PlanResult> = Vec::new();
    let mut transit_ft = 0.0;
    for (i, g) in parts.iter().enumerate() {
        let (plans, _) = bcd_plan(g, p, 4.0);
        all_plans.extend(plans);
        if i > 0 {
            // real inter-parcel hops carry climb + circuit + descent overhead
            // well beyond the straight-line distance between sub-polygon edges
            transit_ft += poly_distance(&parts[i - 1], g) + 1000.0;
        }
    }
    let t = combined_time(&all_plans, p, transit_ft);
    if t < best_t {
        best_plans = all_plans;
        best_t = t;
    }
    (best_plans, best_t)
}
