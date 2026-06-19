//! Boustrophedon cell decomposition + greedy merging — full port of the Python
//! `_x_events` / `_decompose` / `bcd_plan` (+ `_union_plan`).
//!
//! Boolean ops (polygon ∩ slab-box, unary unions) use the `geo` crate
//! (pure-Rust, i_overlay-backed). This is the §5 accuracy-risk zone: float
//! boolean ops can wobble vertices and flip borderline merge decisions, so the
//! merge thresholds and comparison directions match the Python exactly and the
//! result is validated against the reference (tests/golden.rs, single_*).

use std::collections::BTreeMap;

use geo::{BooleanOps, Coord, LineString, MultiPolygon as GMP, Polygon as GPoly};

use crate::field::block_cost;
use crate::geom::{rotate_poly, Geom, Polygon, Pt};
use crate::params::{AircraftParams, ACRE_FT2};
use crate::plan::{coverage, plan_single_direction, PlanResult};
use crate::scanline::{best_single_direction, refine};

// ---------- geo interop ----------

fn to_geo(p: &Polygon) -> GPoly<f64> {
    let close = |r: &[Pt]| {
        let mut v: Vec<Coord<f64>> = r.iter().map(|&[x, y]| Coord { x, y }).collect();
        if let Some(&first) = r.first() {
            v.push(Coord { x: first[0], y: first[1] });
        }
        LineString::new(v)
    };
    GPoly::new(close(&p.exterior), p.interiors.iter().map(|h| close(h)).collect())
}

fn geom_to_geo(g: &Geom) -> Vec<GPoly<f64>> {
    g.parts().iter().map(|p| to_geo(p)).collect()
}

fn from_geo_poly(gp: &GPoly<f64>) -> Polygon {
    let ring = |ls: &LineString<f64>| -> Vec<Pt> {
        let mut v: Vec<Pt> = ls.0.iter().map(|c| [c.x, c.y]).collect();
        // drop the closing duplicate to match our unclosed convention
        if v.len() >= 2 && v[0] == v[v.len() - 1] {
            v.pop();
        }
        v
    };
    Polygon {
        exterior: ring(gp.exterior()),
        interiors: gp.interiors().iter().map(ring).collect(),
    }
}

fn mp_to_geom(mp: GMP<f64>) -> Geom {
    let parts: Vec<Polygon> = mp.0.iter().map(from_geo_poly).collect();
    if parts.len() == 1 {
        Geom::Poly(parts.into_iter().next().unwrap())
    } else {
        Geom::Multi(parts)
    }
}

/// Unary union of a set of Geoms → one Geom (possibly MultiPolygon).
fn union_geoms(geoms: &[&Geom]) -> Geom {
    let polys: Vec<GPoly<f64>> = geoms.iter().flat_map(|g| geom_to_geo(g)).collect();
    mp_to_geom(geo::unary_union(polys.iter()))
}

// ---------- geometry helpers ----------

/// `_x_events`: x-coords where a vertical line's connectivity with the polygon
/// can change = weak local x-extrema on any ring (incl. holes), plus the x
/// bounds. Rounded to 0.1 and de-duplicated, ascending.
fn x_events(poly: &Polygon) -> Vec<f64> {
    let mut xs: Vec<f64> = Vec::new();
    for ring in poly.rings() {
        let n = ring.len();
        if n == 0 {
            continue;
        }
        for i in 0..n {
            let x0 = ring[(i + n - 1) % n][0];
            let x1 = ring[i][0];
            let x2 = ring[(i + 1) % n][0];
            if (x1 <= x0 && x1 <= x2) || (x1 >= x0 && x1 >= x2) {
                xs.push(x1);
            }
        }
    }
    let (minx, _, maxx, _) = poly.bounds();
    xs.push(minx);
    xs.push(maxx);
    let mut r: Vec<f64> = xs.iter().map(|x| (x * 10.0).round() / 10.0).collect();
    r.sort_by(|a, b| a.partial_cmp(b).unwrap());
    r.dedup();
    r
}

/// Total length of collinear, overlapping boundary shared by two geometries.
/// Mirrors `g.boundary.intersection(h.boundary).length` for the axis-aligned
/// cut-line and edge-coincident cases produced by slab decomposition.
fn shared_boundary_len(a: &Geom, b: &Geom) -> f64 {
    let mut total = 0.0;
    for pa in a.parts() {
        for pb in b.parts() {
            total += rings_shared_len(pa, pb);
        }
    }
    total
}

fn rings_shared_len(a: &Polygon, b: &Polygon) -> f64 {
    let edges = |p: &Polygon| -> Vec<(Pt, Pt)> {
        let mut e = Vec::new();
        for ring in p.rings() {
            let n = ring.len();
            for i in 0..n {
                e.push((ring[i], ring[(i + 1) % n]));
            }
        }
        e
    };
    let ea = edges(a);
    let eb = edges(b);
    let mut total = 0.0;
    for &(a0, a1) in &ea {
        for &(b0, b1) in &eb {
            total += collinear_overlap(a0, a1, b0, b1);
        }
    }
    total
}

/// Overlap length of two segments if (near-)collinear and overlapping, else 0.
fn collinear_overlap(a0: Pt, a1: Pt, b0: Pt, b1: Pt) -> f64 {
    let dax = a1[0] - a0[0];
    let day = a1[1] - a0[1];
    let la = dax.hypot(day);
    if la < 1e-9 {
        return 0.0;
    }
    let (ux, uy) = (dax / la, day / la);
    // b endpoints must lie on line a (perpendicular distance ~0)
    let perp = |p: Pt| ((p[0] - a0[0]) * (-uy) + (p[1] - a0[1]) * ux).abs();
    if perp(b0) > 1e-6 || perp(b1) > 1e-6 {
        return 0.0;
    }
    // direction of b must be parallel (collinear already ensured)
    let proj = |p: Pt| (p[0] - a0[0]) * ux + (p[1] - a0[1]) * uy;
    let (mut ta0, mut ta1) = (0.0, la);
    if ta0 > ta1 {
        std::mem::swap(&mut ta0, &mut ta1);
    }
    let (mut tb0, mut tb1) = (proj(b0), proj(b1));
    if tb0 > tb1 {
        std::mem::swap(&mut tb0, &mut tb1);
    }
    let lo = ta0.max(tb0);
    let hi = ta1.min(tb1);
    (hi - lo).max(0.0)
}

/// Min distance between two geometries (part-wise; 0 if any parts touch).
fn geom_distance(a: &Geom, b: &Geom) -> f64 {
    let mut best = f64::INFINITY;
    for pa in a.parts() {
        for pb in b.parts() {
            best = best.min(crate::geom::poly_distance(pa, pb));
        }
    }
    best
}

fn geom_area(g: &Geom) -> f64 {
    g.area()
}

/// `_decompose`: boustrophedon cells. Rotate so sweep ∥ x, slab between events,
/// connected components, absorb slivers into the longest-shared-boundary
/// neighbor, rotate back. Returns cells as Geoms (a merged sliver can be multi).
fn decompose(part: &Polygon, sweep_deg: f64, min_cell_ft2: f64) -> Vec<Geom> {
    let c = part.centroid();
    let rot = rotate_poly(part, -sweep_deg, c);
    let (_, miny, _, maxy) = rot.bounds();
    let xs = x_events(&rot);
    let rot_geo = to_geo(&rot);

    let mut cells: Vec<Geom> = Vec::new();
    for w in xs.windows(2) {
        let (xa, xb) = (w[0], w[1]);
        if xb - xa < 1.0 {
            continue;
        }
        let slab = GPoly::new(
            LineString::new(vec![
                Coord { x: xa, y: miny - 10.0 },
                Coord { x: xb, y: miny - 10.0 },
                Coord { x: xb, y: maxy + 10.0 },
                Coord { x: xa, y: maxy + 10.0 },
                Coord { x: xa, y: miny - 10.0 },
            ]),
            vec![],
        );
        let clipped = rot_geo.intersection(&slab);
        for gp in &clipped.0 {
            let poly = from_geo_poly(gp);
            if poly.area() > 1.0 {
                cells.push(Geom::Poly(poly));
            }
        }
    }

    // absorb slivers into the best-touching neighbor
    let mut changed = true;
    while changed {
        changed = false;
        for i in 0..cells.len() {
            if geom_area(&cells[i]) >= min_cell_ft2 {
                continue;
            }
            let mut best_j: Option<usize> = None;
            let mut best_len = 0.0;
            for j in 0..cells.len() {
                if i == j {
                    continue;
                }
                let l = shared_boundary_len(&cells[i], &cells[j]);
                if l > best_len {
                    best_len = l;
                    best_j = Some(j);
                }
            }
            if let Some(j) = best_j {
                if best_len > 0.5 {
                    let merged = union_geoms(&[&cells[j], &cells[i]]);
                    cells[j] = merged;
                    cells.remove(i);
                    changed = true;
                    break;
                }
            }
        }
    }

    cells.iter().map(|g| rotate_geom(g, sweep_deg, c)).collect()
}

fn rotate_geom(g: &Geom, deg: f64, o: Pt) -> Geom {
    match g {
        Geom::Poly(p) => Geom::Poly(rotate_poly(p, deg, o)),
        Geom::Multi(v) => Geom::Multi(v.iter().map(|p| rotate_poly(p, deg, o)).collect()),
    }
}

/// `_union_plan`: best plan for a merged block, evaluated only at candidate
/// angles, scored by block_cost with the 0.95 coverage gate.
fn union_plan(u: &Geom, angles: &[f64], p: &AircraftParams) -> Option<PlanResult> {
    let area = u.area();
    let mut best: Option<PlanResult> = None;
    for &a in angles {
        let r = plan_single_direction(u, a.rem_euclid(180.0), p);
        if coverage(&r, area, p) >= 0.95 {
            match &best {
                None => best = Some(r),
                Some(b) if block_cost(&r, p) < block_cost(b, p) => best = Some(r),
                _ => {}
            }
        }
    }
    best
}

fn dedup_angles(mut v: Vec<f64>) -> Vec<f64> {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v.dedup();
    v
}

/// `bcd_plan`: decompose at several sweeps, optimal angle per cell, greedy
/// merging, keep the cheapest plan set. Returns (plans, block_cost_total).
pub fn bcd_plan(part: &Polygon, p: &AircraftParams, coarse_deg: f64) -> (Vec<PlanResult>, f64) {
    let pg: Geom = part.clone().into();
    let base_angle = best_single_direction(&pg, p, coarse_deg, 8).angle_deg;
    let base = refine(&pg, base_angle, p);
    let mut best_cost = block_cost(&base, p);
    let mut best_plans = vec![base.clone()];

    let sweeps = dedup_angles(vec![0.0, 90.0, (base.angle_deg.rem_euclid(180.0) * 10.0).round() / 10.0]);
    let min_cell = 0.5 * ACRE_FT2;
    let gap_thresh = p.gap_threshold_ft();

    for sw in sweeps {
        let cells = decompose(part, sw, min_cell);
        if cells.len() <= 1 {
            continue;
        }
        // geoms / plans keyed by id (ordered map for deterministic iteration)
        let mut geoms: BTreeMap<usize, Geom> = BTreeMap::new();
        let mut plans: BTreeMap<usize, PlanResult> = BTreeMap::new();
        for (i, g) in cells.into_iter().enumerate() {
            plans.insert(i, best_single_direction(&g, p, coarse_deg, 8));
            geoms.insert(i, g);
        }
        let mut next_id = geoms.keys().max().copied().unwrap_or(0) + 1;

        let mergeable = |geoms: &BTreeMap<usize, Geom>, i: usize, j: usize| -> bool {
            shared_boundary_len(&geoms[&i], &geoms[&j]) > 0.5
                || geom_distance(&geoms[&i], &geoms[&j]) <= gap_thresh
        };
        let pair_gain = |geoms: &BTreeMap<usize, Geom>,
                         plans: &BTreeMap<usize, PlanResult>,
                         i: usize,
                         j: usize|
         -> Option<(f64, Geom, PlanResult)> {
            let u = union_geoms(&[&geoms[&i], &geoms[&j]]);
            let angs = dedup_angles(vec![plans[&i].angle_deg, plans[&j].angle_deg, 0.0, 90.0]);
            let pu = union_plan(&u, &angs, p)?;
            let gain = block_cost(&plans[&i], p) + block_cost(&plans[&j], p) + p.t_turn_s
                - block_cost(&pu, p);
            Some((gain, u, pu))
        };

        // phase 1: greedy gain-cached merging
        let mut cache: BTreeMap<(usize, usize), Option<(f64, Geom, PlanResult)>> = BTreeMap::new();
        let ids: Vec<usize> = geoms.keys().copied().collect();
        for a in 0..ids.len() {
            for b in (a + 1)..ids.len() {
                let (i, j) = (ids[a], ids[b]);
                if mergeable(&geoms, i, j) {
                    cache.insert((i, j), pair_gain(&geoms, &plans, i, j));
                }
            }
        }
        while geoms.len() > 1 {
            // pick max positive gain (first on tie, by ascending key order)
            let mut best_kv: Option<((usize, usize), f64)> = None;
            for (k, v) in &cache {
                if let Some((gain, _, _)) = v {
                    if *gain > 0.0 && best_kv.map(|(_, g)| *gain > g).unwrap_or(true) {
                        best_kv = Some((*k, *gain));
                    }
                }
            }
            let Some(((i, j), _)) = best_kv else { break };
            let (_, u, pu) = cache[&(i, j)].clone().unwrap();
            let k = next_id;
            next_id += 1;
            geoms.remove(&i);
            geoms.remove(&j);
            plans.remove(&i);
            plans.remove(&j);
            geoms.insert(k, u);
            plans.insert(k, pu);
            cache.retain(|pr, _| pr.0 != i && pr.1 != i && pr.0 != j && pr.1 != j);
            let others: Vec<usize> = geoms.keys().copied().filter(|&m| m != k).collect();
            for m in others {
                if mergeable(&geoms, m, k) {
                    let key = if m < k { (m, k) } else { (k, m) };
                    cache.insert(key, pair_gain(&geoms, &plans, m, k));
                }
            }
        }

        // phase 2: few blocks remain -> full angle search per candidate union
        let mut improved = true;
        while improved && geoms.len() > 1 {
            improved = false;
            let mut best_mv: Option<(f64, usize, usize, Geom, PlanResult)> = None;
            let ids: Vec<usize> = geoms.keys().copied().collect();
            for a in 0..ids.len() {
                for b in (a + 1)..ids.len() {
                    let (i, j) = (ids[a], ids[b]);
                    if !mergeable(&geoms, i, j) {
                        continue;
                    }
                    let u = union_geoms(&[&geoms[&i], &geoms[&j]]);
                    let pu = best_single_direction(&u, p, coarse_deg, 8);
                    let gain = block_cost(&plans[&i], p) + block_cost(&plans[&j], p) + p.t_turn_s
                        - block_cost(&pu, p);
                    if gain > 0.0 && best_mv.as_ref().map(|m| gain > m.0).unwrap_or(true) {
                        best_mv = Some((gain, i, j, u, pu));
                    }
                }
            }
            if let Some((_, i, j, u, pu)) = best_mv {
                let k = next_id;
                next_id += 1;
                geoms.remove(&i);
                geoms.remove(&j);
                plans.remove(&i);
                plans.remove(&j);
                geoms.insert(k, u);
                plans.insert(k, pu);
                improved = true;
            }
        }

        // phase 3: union same-angle (<=5deg) touching/fly-through blocks
        // unconditionally (a real pilot flies them on one continuous grid).
        let ids: Vec<usize> = geoms.keys().copied().collect();
        let mut parent: BTreeMap<usize, usize> = ids.iter().map(|&i| (i, i)).collect();
        fn find(parent: &mut BTreeMap<usize, usize>, x: usize) -> usize {
            let mut x = x;
            while parent[&x] != x {
                let gp = parent[&parent[&x]];
                parent.insert(x, gp);
                x = gp;
            }
            x
        }
        for a in 0..ids.len() {
            for b in (a + 1)..ids.len() {
                let (i, j) = (ids[a], ids[b]);
                let d = (plans[&i].angle_deg - plans[&j].angle_deg).rem_euclid(180.0);
                if d.min(180.0 - d) <= 5.0 && mergeable(&geoms, i, j) {
                    let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                    parent.insert(ri, rj);
                }
            }
        }
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for &i in &ids {
            let r = find(&mut parent, i);
            groups.entry(r).or_default().push(i);
        }
        let mut final_plans: Vec<PlanResult> = Vec::new();
        for members in groups.values() {
            let geom_refs: Vec<&Geom> = members.iter().map(|m| &geoms[m]).collect();
            let u = union_geoms(&geom_refs);
            let lead = *members
                .iter()
                .max_by(|x, y| geom_area(&geoms[x]).partial_cmp(&geom_area(&geoms[y])).unwrap())
                .unwrap();
            final_plans.push(refine(&u, plans[&lead].angle_deg, p));
        }
        let cost: f64 = final_plans.iter().map(|r| block_cost(r, p)).sum::<f64>()
            + (final_plans.len() as f64 - 1.0) * p.t_turn_s;
        if cost < best_cost {
            best_cost = cost;
            best_plans = final_plans;
        }
    }
    (best_plans, best_cost)
}
