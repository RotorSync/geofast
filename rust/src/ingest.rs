//! Ingestion: KML parsing (std, mirrors the Python regex) and lon/lat → feet.
//!
//! The KML parser is dependency-free and faithful to `load_kml`. The projection
//! is a PLACEHOLDER spherical aeqd — it MUST be replaced with an ellipsoidal
//! aeqd (geographiclib-rs) before trusting absolute coordinates; see below.

use crate::geom::{Geom, Polygon, Pt};

/// Parse Polygons (with holes) from a KML file → lon/lat Geom.
/// Mirrors `load_kml`: <Polygon> blocks, outer/inner boundary <coordinates>,
/// tokens split on whitespace then ',' → (lon, lat).
pub fn load_kml(path: &str) -> std::io::Result<Geom> {
    let txt = std::fs::read_to_string(path)?;
    let mut polys: Vec<Polygon> = Vec::new();
    for poly_block in find_blocks(&txt, "<Polygon>", "</Polygon>") {
        let outer = find_blocks(poly_block, "<outerBoundaryIs>", "</outerBoundaryIs>")
            .into_iter()
            .next();
        let Some(outer) = outer else { continue };
        let exterior = parse_ring(outer);
        let interiors: Vec<Vec<Pt>> =
            find_blocks(poly_block, "<innerBoundaryIs>", "</innerBoundaryIs>")
                .into_iter()
                .map(parse_ring)
                .collect();
        polys.push(close_to_open(exterior, interiors));
    }
    Ok(if polys.len() > 1 {
        Geom::Multi(polys)
    } else {
        Geom::Poly(polys.into_iter().next().unwrap_or(Polygon {
            exterior: Vec::new(),
            interiors: Vec::new(),
        }))
    })
}

/// KML rings repeat the first point as last; code_v2.py keeps them closed when
/// building the Shapely Polygon (Shapely auto-closes), but our internal rings are
/// stored unclosed (coords[:-1] convention). Drop a trailing duplicate.
fn close_to_open(mut exterior: Vec<Pt>, mut interiors: Vec<Vec<Pt>>) -> Polygon {
    let drop_dup = |r: &mut Vec<Pt>| {
        if r.len() >= 2 && r[0] == r[r.len() - 1] {
            r.pop();
        }
    };
    drop_dup(&mut exterior);
    for h in &mut interiors {
        drop_dup(h);
    }
    Polygon { exterior, interiors }
}

fn parse_ring(block: &str) -> Vec<Pt> {
    let coords = find_blocks(block, "<coordinates>", "</coordinates>")
        .into_iter()
        .next()
        .unwrap_or("");
    coords
        .split_whitespace()
        .filter_map(|tok| {
            let mut it = tok.split(',');
            let lon = it.next()?.parse::<f64>().ok()?;
            let lat = it.next()?.parse::<f64>().ok()?;
            Some([lon, lat])
        })
        .collect()
}

/// Return the inner text of each `open..close` block (non-overlapping, in order).
fn find_blocks<'a>(txt: &'a str, open: &str, close: &str) -> Vec<&'a str> {
    let mut out = Vec::new();
    let mut rest = txt;
    while let Some(i) = rest.find(open) {
        let after = &rest[i + open.len()..];
        if let Some(j) = after.find(close) {
            out.push(&after[..j]);
            rest = &after[j + close.len()..];
        } else {
            break;
        }
    }
    out
}

const M_PER_FT: f64 = 0.3048; // international foot (PROJ +units=ft)

/// lon/lat → local feet, ellipsoidal azimuthal-equidistant about the geometry
/// centroid. Reproduces pyproj's `+proj=aeqd +units=ft +datum=WGS84` exactly
/// (verified to sub-0.0001 ft): PROJ's ellipsoidal aeqd is the geodesic
/// projection, so we use Karney's inverse geodesic from the center to each
/// point — distance `s` (m) and forward azimuth `α` — then
///   x = s·sin(α)/0.3048,  y = s·cos(α)/0.3048   (x east, y north, in feet).
/// Pure Rust (no PROJ/GEOS), so it cross-compiles to iOS.
/// Returns (geom_ft, center_lonlat) — keep the center for the inverse transform.
pub fn to_local_feet(geom_ll: &Geom) -> (Geom, Pt) {
    use geographiclib_rs::{Geodesic, InverseGeodesic};
    let g = Geodesic::wgs84();
    let c = geom_ll.centroid(); // [lon, lat]
    let (lon0, lat0) = (c[0], c[1]);
    let project = |p: Pt| -> Pt {
        let (lon, lat) = (p[0], p[1]);
        // 4-tuple = (s12 metres, azi1, azi2, a12); azi1 is the forward azimuth.
        let (s12, azi1, _azi2, _a12): (f64, f64, f64, f64) = g.inverse(lat0, lon0, lat, lon);
        let a = azi1.to_radians();
        [s12 * a.sin() / M_PER_FT, s12 * a.cos() / M_PER_FT]
    };
    let map_poly = |poly: &Polygon| Polygon {
        exterior: poly.exterior.iter().map(|&p| project(p)).collect(),
        interiors: poly
            .interiors
            .iter()
            .map(|h| h.iter().map(|&p| project(p)).collect())
            .collect(),
    };
    let out = match geom_ll {
        Geom::Poly(p) => Geom::Poly(map_poly(p)),
        Geom::Multi(v) => Geom::Multi(v.iter().map(map_poly).collect()),
    };
    (out, c)
}

/// Forward projection of a single lon/lat point into the local-feet frame about
/// a given center (the center returned by `to_local_feet`). Mirrors the closure
/// inside `to_local_feet` so a trailer point can be projected into the same
/// frame as the field. Karney inverse geodesic.
pub fn lonlat_to_local_feet(pt: Pt, center_lonlat: Pt) -> Pt {
    use geographiclib_rs::{Geodesic, InverseGeodesic};
    let g = Geodesic::wgs84();
    let (lon0, lat0) = (center_lonlat[0], center_lonlat[1]);
    let (s12, azi1, _azi2, _a12): (f64, f64, f64, f64) = g.inverse(lat0, lon0, pt[1], pt[0]);
    let a = azi1.to_radians();
    [s12 * a.sin() / M_PER_FT, s12 * a.cos() / M_PER_FT]
}

/// Inverse of `to_local_feet` for a single point: local feet → (lon, lat), given
/// the projection center returned by `to_local_feet`. Karney direct geodesic.
pub fn local_feet_to_lonlat(pt: Pt, center_lonlat: Pt) -> Pt {
    use geographiclib_rs::{DirectGeodesic, Geodesic};
    let g = Geodesic::wgs84();
    let s12 = pt[0].hypot(pt[1]) * M_PER_FT;
    let azi1 = pt[0].atan2(pt[1]).to_degrees(); // azimuth from north, east-positive
    let (lat2, lon2): (f64, f64) = g.direct(center_lonlat[1], center_lonlat[0], azi1, s12);
    [lon2, lat2]
}
