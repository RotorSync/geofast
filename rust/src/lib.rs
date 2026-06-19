//! R44 aerial-spray time/cost estimator — Rust port of code_v2.py.
//!
//! Status (see PORTING.md):
//!   DONE  params, geometry primitives, scanline hot path, angle search,
//!         plan bookkeeping, flight simulation, block ordering, KML parse.
//!   STUB  bcd_plan / _decompose (needs i_overlay), ellipsoidal projection
//!         (needs geographiclib-rs), topology-preserving simplify.

mod bind;
pub mod decompose;
pub mod field;
pub mod geom;
pub mod ingest;
pub mod params;
pub mod plan;
pub mod scanline;
pub mod simulate;

use geom::Geom;
use params::AircraftParams;
use plan::PlanResult;
use simulate::{simulate_job, survey_groups, SimResult};

/// `quote`: plan + simulate. Returns (plans, SimResult, trailer_pt).
pub fn quote(
    geom: &Geom,
    p: &AircraftParams,
    trailer: Option<geom::Pt>,
    step_deg: f64,
) -> (Vec<PlanResult>, SimResult, geom::Pt) {
    let (plans, _) = field::plan_field(geom, p, step_deg);
    let n_surveys = survey_groups(
        &geom.parts().into_iter().cloned().collect::<Vec<_>>(),
        p.survey_link_ft,
    );
    let (sim, trailer) = simulate_job(&plans, p, trailer, n_surveys);
    (plans, sim, trailer)
}
