//! Lightweight geometry: dependency-free types + the primitives the hot path
//! needs (rotate, bounds, area, centroid, distances). The cold paths
//! (_decompose / bcd_plan unions) will need a real boolean-op library
//! (i_overlay); see PORTING.md §5.

pub type Pt = [f64; 2];

#[derive(Debug, Clone)]
pub struct Polygon {
    /// Exterior ring, NOT closed (no repeated first==last point), matching the
    /// Python `coords[:-1]` convention used throughout code_v2.py.
    pub exterior: Vec<Pt>,
    pub interiors: Vec<Vec<Pt>>,
}

#[derive(Debug, Clone)]
pub enum Geom {
    Poly(Polygon),
    Multi(Vec<Polygon>),
}

impl Geom {
    pub fn parts(&self) -> Vec<&Polygon> {
        match self {
            Geom::Poly(p) => vec![p],
            Geom::Multi(v) => v.iter().collect(),
        }
    }
    pub fn is_multi(&self) -> bool {
        matches!(self, Geom::Multi(_))
    }

    /// Total filled area over all parts (holes subtracted per part).
    pub fn area(&self) -> f64 {
        self.parts().iter().map(|p| p.area()).sum()
    }

    /// All rings (every part's exterior + interiors). The scanline even-odd rule
    /// handles disjoint parts and holes uniformly (see PORTING.md §3).
    pub fn all_rings(&self) -> Vec<&Vec<Pt>> {
        let mut out = Vec::new();
        for p in self.parts() {
            for r in p.rings() {
                out.push(r);
            }
        }
        out
    }

    /// Concatenated exterior vertices (the `_poly_xy` proxy input). For a
    /// MultiPolygon this is every part's exterior, matching the Python vstack.
    pub fn exterior_verts(&self) -> Vec<Pt> {
        let mut v = Vec::new();
        for p in self.parts() {
            v.extend_from_slice(&p.exterior);
        }
        v
    }

    /// Area-weighted centroid over parts (rotation origin; exact value does not
    /// affect span results — only conditioning — but we match Shapely closely).
    pub fn centroid(&self) -> Pt {
        let parts = self.parts();
        if parts.len() == 1 {
            return parts[0].centroid();
        }
        let mut a_tot = 0.0;
        let mut cx = 0.0;
        let mut cy = 0.0;
        for p in &parts {
            let a = p.area();
            let c = p.centroid();
            a_tot += a;
            cx += c[0] * a;
            cy += c[1] * a;
        }
        if a_tot.abs() < 1e-12 {
            return parts[0].centroid();
        }
        [cx / a_tot, cy / a_tot]
    }
}

impl From<Polygon> for Geom {
    fn from(p: Polygon) -> Self {
        Geom::Poly(p)
    }
}

#[inline]
pub fn dist(a: Pt, b: Pt) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1])
}

impl Polygon {
    /// All rings (exterior first, then interiors), each as an unclosed Vec<Pt>.
    pub fn rings(&self) -> impl Iterator<Item = &Vec<Pt>> {
        std::iter::once(&self.exterior).chain(self.interiors.iter())
    }

    /// Shoelace area of the exterior minus interiors (always positive).
    pub fn area(&self) -> f64 {
        let ring_area = |r: &Vec<Pt>| -> f64 {
            let n = r.len();
            if n < 3 {
                return 0.0;
            }
            let mut s = 0.0;
            for i in 0..n {
                let a = r[i];
                let b = r[(i + 1) % n];
                s += a[0] * b[1] - b[0] * a[1];
            }
            s.abs() / 2.0
        };
        let mut a = ring_area(&self.exterior);
        for hole in &self.interiors {
            a -= ring_area(hole);
        }
        a
    }

    /// (minx, miny, maxx, maxy) over the exterior.
    pub fn bounds(&self) -> (f64, f64, f64, f64) {
        let mut minx = f64::INFINITY;
        let mut miny = f64::INFINITY;
        let mut maxx = f64::NEG_INFINITY;
        let mut maxy = f64::NEG_INFINITY;
        for &[x, y] in &self.exterior {
            minx = minx.min(x);
            miny = miny.min(y);
            maxx = maxx.max(x);
            maxy = maxy.max(y);
        }
        (minx, miny, maxx, maxy)
    }

    /// GEOS-faithful area centroid: exterior moment minus hole moments, divided
    /// by net area (holes contribute negative area). Matches Shapely's
    /// `Polygon.centroid`, which code_v2.py uses as the rotation origin — getting
    /// this bit-close matters because the rotate→clip→rotate-back round-trip's
    /// float rounding depends on the origin, and that can flip near-tie block
    /// entry-config selection in the sim (see PORTING.md §4).
    pub fn centroid(&self) -> Pt {
        let (ecx, ecy, ea) = ring_centroid(&self.exterior);
        let mut mx = ecx * ea;
        let mut my = ecy * ea;
        let mut net = ea;
        for hole in &self.interiors {
            let (hcx, hcy, ha) = ring_centroid(hole);
            mx -= hcx * ha;
            my -= hcy * ha;
            net -= ha;
        }
        if net.abs() < 1e-12 {
            let r = &self.exterior;
            let n = r.len().max(1);
            let (sx, sy) = r.iter().fold((0.0, 0.0), |(sx, sy), p| (sx + p[0], sy + p[1]));
            return [sx / n as f64, sy / n as f64];
        }
        [mx / net, my / net]
    }
}

/// (cx, cy, unsigned_area) of a ring via the shoelace centroid formula.
fn ring_centroid(r: &[Pt]) -> (f64, f64, f64) {
    let n = r.len();
    if n < 3 {
        let (sx, sy) = r.iter().fold((0.0, 0.0), |(sx, sy), p| (sx + p[0], sy + p[1]));
        let d = n.max(1) as f64;
        return (sx / d, sy / d, 0.0);
    }
    let mut a = 0.0;
    let mut cx = 0.0;
    let mut cy = 0.0;
    for i in 0..n {
        let p0 = r[i];
        let p1 = r[(i + 1) % n];
        let cross = p0[0] * p1[1] - p1[0] * p0[1];
        a += cross;
        cx += (p0[0] + p1[0]) * cross;
        cy += (p0[1] + p1[1]) * cross;
    }
    a *= 0.5;
    if a.abs() < 1e-12 {
        let (sx, sy) = r.iter().fold((0.0, 0.0), |(sx, sy), p| (sx + p[0], sy + p[1]));
        return (sx / n as f64, sy / n as f64, 0.0);
    }
    (cx / (6.0 * a), cy / (6.0 * a), a.abs())
}

/// Rotate a single point about origin `o` by `deg` degrees CCW.
///
/// Replicates shapely.affinity.rotate's exact float arithmetic so the
/// rotate→clip→rotate-back round-trip rounds identically to the Python
/// reference: radians via `(deg * PI) / 180` (Python's `math.radians` order, not
/// Rust's `to_radians`), and the affine offset formulation
/// `new = cos*x - sin*y + xoff` rather than `o + cos*dx - sin*dy`. Algebraically
/// equal, but bit-for-bit closer — which matters for near-tie block-config
/// selection (PORTING.md §4).
#[inline]
pub fn rotate_pt(p: Pt, deg: f64, o: Pt) -> Pt {
    let th = (deg * std::f64::consts::PI) / 180.0;
    let cosp = th.cos();
    let sinp = th.sin();
    let xoff = o[0] - o[0] * cosp + o[1] * sinp;
    let yoff = o[1] - o[0] * sinp - o[1] * cosp;
    [cosp * p[0] - sinp * p[1] + xoff, sinp * p[0] + cosp * p[1] + yoff]
}

/// Rotate every ring of a polygon (used to map a planned segment back to world
/// coords, mirroring `rotate(LineString, angle_deg, origin=c)`).
pub fn rotate_poly(poly: &Polygon, deg: f64, o: Pt) -> Polygon {
    let map = |r: &Vec<Pt>| r.iter().map(|&p| rotate_pt(p, deg, o)).collect();
    Polygon {
        exterior: map(&poly.exterior),
        interiors: poly.interiors.iter().map(map).collect(),
    }
}

/// Minimum distance between two polygons. Placeholder: boundary point/segment
/// distance ignoring containment (sufficient for separated parcels, which is the
/// only case code_v2.py uses it on — inter-part ferry & survey linkage). For a
/// fully faithful `Polygon.distance` (which is 0 when one contains the other),
/// add a point-in-polygon test; see PORTING.md §5.
pub fn poly_distance(a: &Polygon, b: &Polygon) -> f64 {
    let mut best = f64::INFINITY;
    for ra in a.rings() {
        for rb in b.rings() {
            best = best.min(ring_ring_dist(ra, rb));
        }
    }
    best
}

fn ring_ring_dist(a: &[Pt], b: &[Pt]) -> f64 {
    let mut best = f64::INFINITY;
    let na = a.len();
    let nb = b.len();
    for i in 0..na {
        let s0 = a[i];
        let s1 = a[(i + 1) % na];
        for j in 0..nb {
            let t0 = b[j];
            let t1 = b[(j + 1) % nb];
            best = best.min(seg_seg_dist(s0, s1, t0, t1));
        }
    }
    best
}

fn seg_seg_dist(p1: Pt, p2: Pt, p3: Pt, p4: Pt) -> f64 {
    // If they intersect, distance is 0; otherwise min of endpoint-to-segment.
    if segments_intersect(p1, p2, p3, p4) {
        return 0.0;
    }
    let d = [
        pt_seg_dist(p1, p3, p4),
        pt_seg_dist(p2, p3, p4),
        pt_seg_dist(p3, p1, p2),
        pt_seg_dist(p4, p1, p2),
    ];
    d.iter().cloned().fold(f64::INFINITY, f64::min)
}

fn pt_seg_dist(p: Pt, a: Pt, b: Pt) -> f64 {
    let vx = b[0] - a[0];
    let vy = b[1] - a[1];
    let wx = p[0] - a[0];
    let wy = p[1] - a[1];
    let len2 = vx * vx + vy * vy;
    let t = if len2 > 0.0 {
        ((wx * vx + wy * vy) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let cx = a[0] + t * vx;
    let cy = a[1] + t * vy;
    (p[0] - cx).hypot(p[1] - cy)
}

fn segments_intersect(p1: Pt, p2: Pt, p3: Pt, p4: Pt) -> bool {
    let d = |a: Pt, b: Pt, c: Pt| (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    let d1 = d(p3, p4, p1);
    let d2 = d(p3, p4, p2);
    let d3 = d(p1, p2, p3);
    let d4 = d(p1, p2, p4);
    ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0))
        && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0))
}
