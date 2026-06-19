# geofast-native

The spray line-generation engine for geofast, built as a Python extension module
(`geofast._native`) via [PyO3](https://pyo3.rs) + [maturin](https://www.maturin.rs).

This is the RotorSync **field-predictor** planner: given a field boundary it does
an ellipsoidal azimuthal-equidistant projection, a sweep-angle search, per-gap
fly-through-vs-turn decisions, and boustrophedon cell decomposition to produce the
spray pass lines.

## Layout

- `src/geom.rs`, `src/params.rs` — geometry primitives and aircraft/swath constants
- `src/ingest.rs` — lon/lat ⇄ local-feet projection (pure-Rust geodesic, no PROJ/GEOS)
- `src/scanline.rs`, `src/plan.rs`, `src/field.rs`, `src/decompose.rs`, `src/simulate.rs` — the planner
- `src/lib.rs` — crate root + `quote()` entry point
- `src/bind.rs` — the PyO3 binding (`plan_lines`) called by `geofast/spray_line_generator.py`

The dependency stack (`geo`, `geographiclib-rs`, transitively `i_overlay`) is pure
Rust with no C dependencies, so it builds into a self-contained wheel.

## Build

Built automatically by maturin from the repo root (`pyproject.toml` points its
`manifest-path` here). For local development:

```sh
maturin develop --release   # compile + install the extension into the active venv
```

End users installing a prebuilt wheel do **not** need a Rust toolchain; only
building from source does.
