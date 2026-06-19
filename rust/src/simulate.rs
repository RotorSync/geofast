//! Flight simulation: block ordering, entry/exit configs, position-true reloads.
//! Pure-logic port of the Python (no geometry library needed). The tie-break
//! parity notes in PORTING.md §4 apply to `_order_blocks` and the trailer default.

use crate::geom::{dist, Pt};
use crate::params::{AircraftParams, ACRE_FT2};
use crate::plan::{LinePart, PlanResult};

#[derive(Debug, Clone, Default)]
pub struct SimResult {
    pub setup_s: f64,
    pub spray_s: f64,
    pub turn_s: f64,
    pub dead_s: f64,
    pub ferry_s: f64,
    pub load_s: f64,
    pub n_loads: i64,
    pub n_turns: i64,
    pub dead_ft: f64,
    pub ferry_ft: f64,
}

impl SimResult {
    pub fn total_s(&self) -> f64 {
        self.setup_s + self.spray_s + self.turn_s + self.dead_s + self.ferry_s + self.load_s
    }
}

/// One traversed part with its source line index: ((p0, p1, is_spray, len), li).
type Step = (LinePart, usize);

/// `_traverse`: serpentine sequence for a block under entry config
/// (reverse_lines, start_back).
pub fn traverse(plan: &PlanResult, reverse_lines: bool, start_back: bool) -> Vec<Step> {
    let mut seq = Vec::new();
    let n = plan.lines.len();
    let mut backwards = start_back;
    for idx in 0..n {
        let li = if reverse_lines { n - 1 - idx } else { idx };
        let parts = &plan.lines[li];
        if backwards {
            for &(a, b, sp, l) in parts.iter().rev() {
                seq.push(((b, a, sp, l), li));
            }
        } else {
            for &(a, b, sp, l) in parts.iter() {
                seq.push(((a, b, sp, l), li));
            }
        }
        backwards = !backwards;
    }
    seq
}

/// `_block_traversal_stats`: (time_excl_loads, entry_pt, exit_pt, n_turns).
fn block_stats(plan: &PlanResult, cfg: (bool, bool), p: &AircraftParams) -> (f64, Option<Pt>, Option<Pt>, i64) {
    let seq = traverse(plan, cfg.0, cfg.1);
    if seq.is_empty() {
        return (0.0, None, None, 0);
    }
    let t_len: f64 = seq.iter().map(|((_, _, _, l), _)| *l).sum::<f64>() / p.speed_fps;
    // n_turns = number of distinct line indices - 1
    let mut lis: Vec<usize> = seq.iter().map(|(_, li)| *li).collect();
    lis.sort_unstable();
    lis.dedup();
    let n_turns = lis.len() as i64 - 1;
    let t = t_len + n_turns as f64 * p.t_turn_s;
    let entry = (seq[0].0).0;
    let exit = (seq[seq.len() - 1].0).1;
    (t, Some(entry), Some(exit), n_turns)
}

const CONFIGS: [(bool, bool); 4] = [(false, false), (false, true), (true, false), (true, true)];

/// Lexicographic permutations of 0..k (matches itertools.permutations order),
/// via repeated next_permutation starting from the sorted sequence.
fn permutations(k: usize) -> Vec<Vec<usize>> {
    let mut cur: Vec<usize> = (0..k).collect();
    let mut out = vec![cur.clone()];
    while next_permutation(&mut cur) {
        out.push(cur.clone());
    }
    out
}

/// In-place next lexicographic permutation; false when already the last.
fn next_permutation(a: &mut [usize]) -> bool {
    let n = a.len();
    if n < 2 {
        return false;
    }
    let mut i = n - 1;
    while i > 0 && a[i - 1] >= a[i] {
        i -= 1;
    }
    if i == 0 {
        return false;
    }
    let mut j = n - 1;
    while a[j] <= a[i - 1] {
        j -= 1;
    }
    a.swap(i - 1, j);
    a[i..].reverse();
    true
}

type BlockStat = (f64, Option<Pt>, Option<Pt>, i64);

/// `_order_blocks`: choose block order + entry configs. Exact (perm + DP over
/// configs) for k<=max_exact, greedy otherwise. Returns ordered (block_index, cfg).
pub fn order_blocks(
    plans: &[PlanResult],
    trailer: Pt,
    p: &AircraftParams,
    max_exact: usize,
) -> Vec<(usize, (bool, bool))> {
    let k = plans.len();
    // stats[block][config]
    let stats: Vec<[BlockStat; 4]> = plans
        .iter()
        .map(|pl| {
            [
                block_stats(pl, CONFIGS[0], p),
                block_stats(pl, CONFIGS[1], p),
                block_stats(pl, CONFIGS[2], p),
                block_stats(pl, CONFIGS[3], p),
            ]
        })
        .collect();

    let entry = |b: usize, c: usize| stats[b][c].1.unwrap_or([0.0, 0.0]);
    let exit = |b: usize, c: usize| stats[b][c].2.unwrap_or([0.0, 0.0]);
    let tcost = |b: usize, c: usize| stats[b][c].0;

    // cost of a full permutation via DP over the 4 entry configs
    let seq_cost = |perm: &[usize]| -> f64 {
        let inf = f64::INFINITY;
        let first = perm[0];
        let mut dp: Vec<f64> = (0..4)
            .map(|ci| dist(trailer, entry(first, ci)) / p.ferry_fps + tcost(first, ci))
            .collect();
        for w in perm.windows(2) {
            let (prev, cur) = (w[0], w[1]);
            let mut row = vec![0.0; 4];
            for ci in 0..4 {
                let mut best = inf;
                for pj in 0..4 {
                    let tr = dist(exit(prev, pj), entry(cur, ci)) / p.speed_fps + p.t_turn_s;
                    best = best.min(dp[pj] + tr);
                }
                row[ci] = best + tcost(cur, ci);
            }
            dp = row;
        }
        let last = perm[perm.len() - 1];
        (0..4)
            .map(|ci| dp[ci] + dist(exit(last, ci), trailer) / p.ferry_fps)
            .fold(inf, f64::min)
    };

    let perms: Vec<Vec<usize>> = if k == 1 {
        vec![vec![0]]
    } else if k <= max_exact {
        permutations(k)
    } else {
        // greedy nearest-entry from the trailer
        let mut left: Vec<usize> = (0..k).collect();
        let mut order = Vec::new();
        let mut pos = trailer;
        while !left.is_empty() {
            // min over remaining of min over configs of dist(pos, entry)
            let (idx_in_left, _) = left
                .iter()
                .enumerate()
                .min_by(|(_, &i), (_, &j)| {
                    let di = (0..4).map(|c| dist(pos, entry(i, c))).fold(f64::INFINITY, f64::min);
                    let dj = (0..4).map(|c| dist(pos, entry(j, c))).fold(f64::INFINITY, f64::min);
                    di.partial_cmp(&dj).unwrap()
                })
                .unwrap();
            let nxt = left.remove(idx_in_left);
            order.push(nxt);
            pos = exit(nxt, 0); // stats[nxt][0][2]
        }
        vec![order]
    };

    // pick first minimizing permutation (first-on-tie parity)
    let mut best_perm = &perms[0];
    let mut best_c = seq_cost(best_perm);
    for perm in perms.iter().skip(1) {
        let c = seq_cost(perm);
        if c < best_c {
            best_c = c;
            best_perm = perm;
        }
    }

    // exact config recovery: DP with backpointers along the winning perm
    let first = best_perm[0];
    let mut dp: Vec<f64> = (0..4)
        .map(|ci| dist(trailer, entry(first, ci)) / p.ferry_fps + tcost(first, ci))
        .collect();
    let mut back: Vec<[usize; 4]> = vec![[0; 4]];
    for w in best_perm.windows(2) {
        let (prev, cur) = (w[0], w[1]);
        let mut row = vec![0.0; 4];
        let mut brow = [0usize; 4];
        for ci in 0..4 {
            let mut best_v = f64::INFINITY;
            let mut best_pj = 0;
            for pj in 0..4 {
                let v = dp[pj] + dist(exit(prev, pj), entry(cur, ci)) / p.speed_fps + p.t_turn_s;
                if v < best_v {
                    best_v = v;
                    best_pj = pj;
                }
            }
            row[ci] = best_v + tcost(cur, ci);
            brow[ci] = best_pj;
        }
        dp = row;
        back.push(brow);
    }
    let last = best_perm[best_perm.len() - 1];
    let mut ci = (0..4)
        .min_by(|&a, &b| {
            let va = dp[a] + dist(exit(last, a), trailer) / p.ferry_fps;
            let vb = dp[b] + dist(exit(last, b), trailer) / p.ferry_fps;
            va.partial_cmp(&vb).unwrap()
        })
        .unwrap();
    let mut cfgs = vec![ci];
    for lvl in (1..best_perm.len()).rev() {
        ci = back[lvl][ci];
        cfgs.push(ci);
    }
    cfgs.reverse();

    best_perm
        .iter()
        .zip(cfgs.iter())
        .map(|(&bi, &c)| (bi, CONFIGS[c]))
        .collect()
}

/// `survey_groups`: single-linkage clustering of parts within link_ft.
pub fn survey_groups(parts: &[crate::geom::Polygon], link_ft: f64) -> i64 {
    let n = parts.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if crate::geom::poly_distance(&parts[i], &parts[j]) <= link_ft {
                let pi = find(&mut parent, i);
                let pj = find(&mut parent, j);
                parent[pi] = pj;
            }
        }
    }
    let mut roots: Vec<usize> = (0..n).map(|i| find(&mut parent, i)).collect();
    roots.sort_unstable();
    roots.dedup();
    roots.len() as i64
}

/// `simulate_job`: position-true flight sim. Returns (SimResult, trailer_pt).
///
/// Block ORDER comes from `order_blocks` (permutation + DP heuristic). The
/// per-block entry CONFIG, however, is chosen by *actually simulating* the
/// candidate configs and taking the lowest total time, for small block counts.
/// The 4 configs are identical in spray/turn/dead but differ in entry/exit and
/// in where the position-true reloads happen; `order_blocks` scores configs on
/// entry/exit ferry only, so near-ties there flip on ~1e-9 scanline-vs-GEOS
/// coordinate noise and shift reload ferry (PORTING.md §4). Brute-forcing the
/// config by real sim removes that sensitivity and is deterministic and optimal
/// — Rust is never worse than the reference. Falls back to the heuristic config
/// when the block count makes 4^k enumeration impractical.
pub fn simulate_job(
    plans: &[PlanResult],
    p: &AircraftParams,
    trailer: Option<Pt>,
    n_surveys: i64,
) -> (SimResult, Pt) {
    let trailer = trailer.unwrap_or_else(|| default_trailer(plans));
    let ordered = order_blocks(plans, trailer, p, 6);
    let setup_s = p.t_setup_s + n_surveys as f64 * p.t_survey_s;

    // Block entry-config selection (PORTING.md §4):
    //  - OPTIMAL_CONFIG = true  (chosen default): brute-force all 4^k configs for
    //    the chosen order by full sim and take the lowest total. Deterministic
    //    and globally optimal — removes the ~1e-9-noise tie-break sensitivity.
    //    Eliminates the single_08-style +1.5% outlier; on a few fields the plan
    //    then beats (under-quotes vs) the calibrated Python reference by ≤1.7%.
    //  - OPTIMAL_CONFIG = false: faithful mode — use order_blocks' order+configs
    //    as-is to reproduce Python exactly, accepting the rare reload-ferry
    //    tie-break flip.
    const OPTIMAL_CONFIG: bool = true;

    let seq: Vec<(usize, (bool, bool))> = if OPTIMAL_CONFIG && !ordered.is_empty() && ordered.len() <= 4 {
        let order: Vec<usize> = ordered.iter().map(|(bi, _)| *bi).collect();
        let k = order.len();
        let mut best: Option<(f64, Vec<(bool, bool)>)> = None;
        for mask in 0..(1usize << (2 * k)) {
            let cfgs: Vec<(bool, bool)> =
                (0..k).map(|i| CONFIGS[(mask >> (2 * i)) & 0b11]).collect();
            let s: Vec<_> = order.iter().copied().zip(cfgs.iter().copied()).collect();
            let t = run_sim(&s, plans, p, trailer, setup_s).total_s();
            if best.as_ref().map(|(bt, _)| t < *bt).unwrap_or(true) {
                best = Some((t, cfgs));
            }
        }
        order.into_iter().zip(best.unwrap().1).collect()
    } else {
        ordered
    };

    let res = run_sim(&seq, plans, p, trailer, setup_s);
    (res, trailer)
}

/// Run the position-true sim for a fixed ordered list of (block, config).
fn run_sim(
    ordered: &[(usize, (bool, bool))],
    plans: &[PlanResult],
    p: &AircraftParams,
    trailer: Pt,
    setup_s: f64,
) -> SimResult {
    let gal_per_ft = p.swath_ft * p.gpa / ACRE_FT2;
    let mut res = SimResult {
        setup_s,
        n_loads: 1, // first load at startup
        ..Default::default()
    };
    let mut tank = p.tank_gal;
    let mut pos = trailer;

    for &(bi, cfg) in ordered {
        let plan = &plans[bi];
        let seq = traverse(plan, cfg.0, cfg.1);
        if seq.is_empty() {
            continue;
        }
        // second turn per line break (traversal charges one at the transition)
        res.turn_s += plan.n_breaks as f64 * p.t_turn_s;
        res.n_turns += plan.n_breaks;
        let entry = (seq[0].0).0;
        let d = dist(pos, entry);
        if pos == trailer {
            res.ferry_ft += d;
            res.ferry_s += d / p.ferry_fps;
        } else {
            res.dead_ft += d;
            res.dead_s += d / p.dead_fps + p.t_turn_s;
            res.n_turns += 1;
        }
        let mut last_li = seq[0].1;
        for ((a, b, is_spray, l), li) in seq {
            if li != last_li {
                // headland turn: clock time + outside-field arc traversal
                res.turn_s += p.t_turn_s + p.arc_ft_per_turn / p.ferry_fps;
                res.n_turns += 1;
                let extra = dist(pos, a) - 2.0 * p.swath_ft;
                if extra > 0.0 {
                    res.dead_ft += extra;
                    res.dead_s += extra / p.dead_fps;
                }
                last_li = li;
            }
            if !is_spray {
                res.dead_ft += l;
                res.dead_s += l / p.dead_fps;
                pos = b;
                continue;
            }
            let mut remaining = l;
            let (ux, uy) = if l > 0.0 {
                ((b[0] - a[0]) / l, (b[1] - a[1]) / l)
            } else {
                (0.0, 0.0)
            };
            let mut cur = a;
            while remaining > 1e-9 {
                let can_ft = if gal_per_ft > 0.0 { tank / gal_per_ft } else { remaining };
                let fly = remaining.min(can_ft);
                res.spray_s += fly / p.speed_fps;
                tank -= fly * gal_per_ft;
                cur = [cur[0] + ux * fly, cur[1] + uy * fly];
                remaining -= fly;
                if remaining > 1e-9 {
                    // reload at cur
                    let dd = dist(cur, trailer);
                    res.ferry_ft += dd;
                    res.ferry_s += dd / p.ferry_fps;
                    res.load_s += p.t_load_s + p.t_cycle_s;
                    res.n_loads += 1;
                    tank = p.tank_gal;
                    let db = dist(trailer, cur);
                    res.ferry_ft += db;
                    res.ferry_s += db / p.ferry_fps;
                }
            }
            pos = b;
        }
    }
    // final return to trailer
    let d = dist(pos, trailer);
    res.ferry_ft += d;
    res.ferry_s += d / p.ferry_fps;
    res
}

/// Default trailer: nearest line-part endpoint to (minx-1, miny-1) over all runs.
/// (Mirrors the Python `unary_union(... buffer(1) ...).bounds[:2]`, whose axis
/// bounds are exactly the global min coords minus the 1 ft buffer.)
fn default_trailer(plans: &[PlanResult]) -> Pt {
    let mut bx = f64::INFINITY;
    let mut by = f64::INFINITY;
    for pl in plans {
        for &(a, b, _) in &pl.runs {
            bx = bx.min(a[0]).min(b[0]);
            by = by.min(a[1]).min(b[1]);
        }
    }
    bx -= 1.0;
    by -= 1.0;
    let target = [bx, by];
    let mut best = [0.0, 0.0];
    let mut best_d = f64::INFINITY;
    for pl in plans {
        for ln in &pl.lines {
            for prt in ln {
                for c in [prt.0, prt.1] {
                    let dd = dist(c, target);
                    if dd < best_d {
                        best_d = dd;
                        best = c;
                    }
                }
            }
        }
    }
    best
}
