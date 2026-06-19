//! PlanResult + plan_single_direction. The PlanResult bookkeeping is a direct
//! port; plan_single_direction consumes the scanline parts from scanline.rs.

use crate::geom::{rotate_pt, Geom, Pt};
use crate::params::{AircraftParams, ACRE_FT2};

/// One per-scanline part after gap classification: a spray or boom-off segment.
pub type LinePart = (Pt, Pt, bool, f64); // (p0, p1, is_spray, length)

#[derive(Debug, Clone, Default)]
pub struct PlanResult {
    pub angle_deg: f64,
    pub spray_ft: f64,
    pub dead_ft: f64,
    pub n_runs: i64,
    /// For plotting only: (p0, p1, is_spray) of every emitted segment.
    pub runs: Vec<(Pt, Pt, bool)>,
    /// Per scanline (after big-gap breaks): list of LinePart. Drives the sim.
    pub lines: Vec<Vec<LinePart>>,
    pub n_breaks: i64,
    pub area_ac: f64,
}

impl PlanResult {
    pub fn time_s(&self, p: &AircraftParams) -> f64 {
        let t_air = self.spray_ft / p.speed_fps + self.dead_ft / p.dead_fps;
        let t_turns = (self.n_runs - 1).max(0) as f64 * p.t_turn_s;
        let n_loads = (self.area_ac * p.gpa / p.tank_gal).ceil();
        let t_loads = n_loads * p.t_load_s;
        let t_ferry = 2.0 * n_loads * p.ferry_ft / p.ferry_fps;
        t_air + t_turns + t_loads + t_ferry
    }
}

/// Mirrors the Python coverage gate `r.spray_ft * p.swath_ft / max(area, 1.0)`.
#[inline]
pub fn coverage(r: &PlanResult, poly_area: f64, p: &AircraftParams) -> f64 {
    r.spray_ft * p.swath_ft / poly_area.max(1.0)
}

/// Build a PlanResult from the per-scanline parts produced by
/// `segments_for_angle`. `parts_per_line` mirrors the Python `lines` list where
/// each entry is `(y, parts)` and `parts` may contain `None` (a break) encoded
/// here as `Part::Break`.
pub fn plan_single_direction(
    geom: &Geom,
    angle_deg: f64,
    p: &AircraftParams,
) -> PlanResult {
    let (lines, c) = crate::scanline::segments_for_angle(geom, angle_deg, p);
    let mut res = PlanResult {
        angle_deg,
        area_ac: geom.area() / ACRE_FT2,
        ..Default::default()
    };
    for (y, parts) in lines {
        let mut run_open = false;
        let mut line_parts: Vec<LinePart> = Vec::new();
        for part in parts {
            match part {
                Part::Break => {
                    run_open = false;
                    if !line_parts.is_empty() {
                        res.lines.push(std::mem::take(&mut line_parts));
                        res.n_breaks += 1;
                    }
                }
                Part::Seg(x0, x1, is_spray) => {
                    // rotate (x0,y)-(x1,y) back to world by +angle about c
                    let a = rotate_pt([x0, y], angle_deg, c);
                    let b = rotate_pt([x1, y], angle_deg, c);
                    res.runs.push((a, b, is_spray));
                    let l = x1 - x0;
                    line_parts.push((a, b, is_spray, l));
                    if is_spray {
                        res.spray_ft += l;
                        if !run_open {
                            res.n_runs += 1;
                            run_open = true;
                        }
                    } else {
                        res.dead_ft += l;
                    }
                }
            }
        }
        if !line_parts.is_empty() {
            res.lines.push(line_parts);
        }
    }
    res
}

/// A classified part on a scanline (in the rotated frame): either a segment
/// `(x0, x1, is_spray)` or a `Break` (too-big gap → turn around / new line).
#[derive(Debug, Clone, Copy)]
pub enum Part {
    Seg(f64, f64, bool),
    Break,
}
