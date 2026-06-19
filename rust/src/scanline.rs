//! The HOT path. Replaces Shapely/GEOS `rot.intersection(scans)` with a direct
//! scanline polygon clip, plus the GEOS-free proxy angle search.
//!
//! Parity contract (see PORTING.md §3): for any angle this must produce the same
//! set of (x0,x1) spans per scanline as the Python GEOS intersection, to ~1e-9.

use crate::geom::{rotate_pt, Geom, Pt};
use crate::params::AircraftParams;
use crate::plan::{coverage, plan_single_direction, Part, PlanResult};

/// Rotate every ring of the geometry into the sweep-aligned frame (Shapely
/// `rotate(geom, -angle_deg, origin=c)`), returning rotated rings (exterior of
/// ring 0 used for bounds) and the centroid used as the rotation origin.
fn rotate_into_frame(geom: &Geom, angle_deg: f64) -> (Vec<Vec<Pt>>, Pt) {
    let c = geom.centroid();
    let rings: Vec<Vec<Pt>> = geom
        .all_rings()
        .iter()
        .map(|r| r.iter().map(|&p| rotate_pt(p, -angle_deg, c)).collect())
        .collect();
    (rings, c)
}

/// All inside-spans of `rings` cut by the horizontal line y = `y`.
/// Even-odd rule over crossings of every ring. Half-open edge convention
/// (`min(y0,y1) <= y < max(y0,y1)`) so a scanline through a vertex is counted
/// once. Returns sorted (x0, x1) pairs; zero-length spans are dropped to mirror
/// `_extract_lines` ignoring degenerate grazes.
fn scanline_spans(rings: &[Vec<Pt>], y: f64) -> Vec<(f64, f64)> {
    let mut xs: Vec<f64> = Vec::new();
    for ring in rings {
        let n = ring.len();
        if n < 2 {
            continue;
        }
        for i in 0..n {
            let a = ring[i];
            let b = ring[(i + 1) % n];
            let (y0, y1) = (a[1], b[1]);
            if y0 == y1 {
                continue; // horizontal edge: no crossing
            }
            // half-open: lower-inclusive, upper-exclusive (orientation-independent)
            let ylo = y0.min(y1);
            let yhi = y0.max(y1);
            if y >= ylo && y < yhi {
                let x = a[0] + (y - y0) * (b[0] - a[0]) / (y1 - y0);
                xs.push(x);
            }
        }
    }
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut spans = Vec::with_capacity(xs.len() / 2);
    let mut i = 0;
    while i + 1 < xs.len() {
        let (x0, x1) = (xs[i], xs[i + 1]);
        if x1 - x0 > 1e-6 {
            spans.push((x0, x1));
        }
        i += 2;
    }
    spans
}

/// Port of `_segments_for_angle`: returns per-scanline `(y, parts)` in the
/// rotated frame plus the centroid `c`. Parts are gap-classified left→right.
pub fn segments_for_angle(
    geom: &Geom,
    angle_deg: f64,
    p: &AircraftParams,
) -> (Vec<(f64, Vec<Part>)>, Pt) {
    let (rings, c) = rotate_into_frame(geom, angle_deg);
    // bounds of the whole rotated geometry
    let mut miny = f64::INFINITY;
    let mut maxy = f64::NEG_INFINITY;
    for ring in &rings {
        for &[_, y] in ring {
            miny = miny.min(y);
            maxy = maxy.max(y);
        }
    }
    // scanline y's: while y <= maxy: y += swath, starting at miny + swath/2
    let mut ys = Vec::new();
    let mut y = miny + p.swath_ft / 2.0;
    while y <= maxy {
        ys.push(y);
        y += p.swath_ft;
    }
    if ys.is_empty() {
        return (Vec::new(), c);
    }

    let thresh = p.gap_threshold_ft();
    let mut lines_out = Vec::with_capacity(ys.len());
    for yy in ys {
        let spans = scanline_spans(&rings, yy);
        if spans.is_empty() {
            continue;
        }
        // walk segments left->right, classify gaps (mirrors the Python walk)
        let mut parts: Vec<Part> = Vec::new();
        let mut prev_end: Option<f64> = None;
        for (x0, x1) in spans {
            if let Some(pe) = prev_end {
                let gap = x0 - pe;
                if gap <= thresh {
                    parts.push(Part::Seg(pe, x0, false)); // fly through, spray off
                } else {
                    parts.push(Part::Break);
                }
            }
            parts.push(Part::Seg(x0, x1, true));
            prev_end = Some(x1);
        }
        lines_out.push((yy, parts));
    }
    (lines_out, c)
}

/// `_angle_scores`: GEOS-free per-angle proxy = spray_time + turn_time, via
/// vertex projection onto each angle's normal. `verts` is the exterior vertex
/// set (use all parts' exteriors for a multipolygon).
pub fn angle_scores(verts: &[Pt], area: f64, angles_deg: &[f64], p: &AircraftParams) -> Vec<f64> {
    let spray_time = (area / p.swath_ft) / p.speed_fps;
    angles_deg
        .iter()
        .map(|&a| {
            let th = a.to_radians();
            let (s, cth) = th.sin_cos();
            // normal n = (-sin, cos); projection of each vertex onto n
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for &[vx, vy] in verts {
                let proj = vx * (-s) + vy * cth;
                lo = lo.min(proj);
                hi = hi.max(proj);
            }
            let width_perp = hi - lo;
            let n_pass = (width_perp / p.swath_ft).ceil();
            let turn_time = (n_pass - 1.0).max(0.0) * p.t_turn_s;
            spray_time + turn_time
        })
        .collect()
}

/// `best_single_direction`: proxy rank → diversified top-K + cardinals → real
/// clip → ±step refine. `step_deg` default 2.0, `top_k` default 8.
pub fn best_single_direction(
    geom: &Geom,
    p: &AircraftParams,
    step_deg: f64,
    top_k: usize,
) -> PlanResult {
    let area = geom.area();
    // angles = arange(0,180,step)
    let mut angles = Vec::new();
    let mut a = 0.0;
    while a < 180.0 {
        angles.push(a);
        a += step_deg;
    }
    let verts = geom.exterior_verts();
    let proxy = angle_scores(&verts, area, &angles, p);

    // STABLE argsort (parity with numpy argsort) — see PORTING.md §4.
    let mut order: Vec<usize> = (0..angles.len()).collect();
    order.sort_by(|&i, &j| proxy[i].partial_cmp(&proxy[j]).unwrap());

    // diversify: skip angles within 8 deg (circular on 180) of a chosen one
    let mut chosen: Vec<usize> = Vec::new();
    for &ai in &order {
        let av = angles[ai];
        let ok = chosen.iter().all(|&c| {
            let d = (av - angles[c]).abs();
            d.min(180.0 - d) >= 8.0
        });
        if ok {
            chosen.push(ai);
        }
        if chosen.len() >= top_k {
            break;
        }
    }
    // always include cardinals 0 and 90
    for card in [0.0_f64, 90.0] {
        let ci = ((card / step_deg).round() as usize) % angles.len();
        if !chosen.contains(&ci) {
            chosen.push(ci);
        }
    }

    let mut best: Option<PlanResult> = None;
    for &ai in &chosen {
        let r = plan_single_direction(geom, angles[ai], p);
        let cov = coverage(&r, area, p);
        if cov >= 0.95 {
            match &best {
                None => best = Some(r),
                Some(b) if r.time_s(p) < b.time_s(p) => best = Some(r),
                _ => {}
            }
        }
    }
    let mut best = match best {
        Some(b) => b,
        None => return plan_single_direction(geom, 0.0, p),
    };
    // local refine ±step around the winner
    for da in [-step_deg, step_deg] {
        let ang = (best.angle_deg + da).rem_euclid(180.0);
        let r = plan_single_direction(geom, ang, p);
        let cov = coverage(&r, area, p);
        if cov >= 0.95 && r.time_s(p) < best.time_s(p) {
            best = r;
        }
    }
    best
}

/// `_refine`: fine 1-degree refinement over [-4, +4] around a coarse best angle.
pub fn refine(geom: &Geom, a0: f64, p: &AircraftParams) -> PlanResult {
    let area = geom.area();
    let mut best: Option<PlanResult> = None;
    for da in -4..=4 {
        let ang = (a0 + da as f64).rem_euclid(180.0);
        let r = plan_single_direction(geom, ang, p);
        let cov = coverage(&r, area, p);
        if cov >= 0.95 {
            match &best {
                None => best = Some(r),
                Some(b) if r.time_s(p) < b.time_s(p) => best = Some(r),
                _ => {}
            }
        }
    }
    best.unwrap_or_else(|| plan_single_direction(geom, a0, p))
}
