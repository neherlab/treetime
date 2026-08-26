# Coalescent output uses a per-segment schema

v1 serializes the inferred coalescent time scale from the `timetree` command with its own per-segment schema. This differs from v0, which writes one interpolated row per date.

**Type**: Output-format change from v0.

**Status**: Approved by the maintainer team.

**v0**: `MergerModel` writes a single `skyline.tsv` with columns `date`, `N_e`, `lower`, `upper`, one row per interpolation date across the tree's span ([`packages/legacy/treetime/treetime/CLI_io.py#L248`](../../packages/legacy/treetime/treetime/CLI_io.py#L248)). The coalescent time scale $T_c$ is not reported; only the derived effective population size $N_e$ and its band appear.

**v1**: The serializers in [`packages/treetime/src/commands/timetree/output/coalescent.rs`](../../packages/treetime/src/commands/timetree/output/coalescent.rs) emit two shapes from one document:

- Rich JSON (`timetree.coalescent.json`): a complete `inputs` block (mode, the skyline `n_points` and `stiffness` where they apply, the confidence width, and generations per year) and an `outputs` block (the coalescent `log_likelihood` where inferred, and `segments`), each segment carrying a numeric-date `segment` interval and nested `T_c` and `N_e` objects with `value`/`lower`/`upper`. The JSON carries the full document; a field that does not apply to the mode (for example `n_points` for a constant or fixed `T_c`, or `log_likelihood` for a fixed `T_c`) is omitted.
- Flat delimited TSV/CSV (`timetree.coalescent.tsv`, `timetree.coalescent.csv`): the selected per-segment projection, one row per segment with dotted columns `segment.start`, `segment.end`, `T_c.value`, `T_c.lower`, `T_c.upper`, `N_e.value`, `N_e.lower`, `N_e.upper`, and no inputs or likelihood. The TSV is in the default `--output-all` set; CSV and JSON are opt-in.

## Decision

Report the coalescent as per-segment intervals rather than per-date points, and report both $T_c$ and the derived $N_e = T_c \cdot \text{gen\_per\_year}$ with their confidence bands. The rich JSON is the single source of truth; the flat rows are its projection, shifting the 0-based JSON index to 1-based.

A fixed, user-supplied $T_c$ writes one band-less segment spanning the tree, carrying no confidence band and no likelihood; the segment span is the range of the coalescent breakpoints. The stderr $N_e$ report uses the same gate, so the file and the screen agree. This diverges from v0, which writes nothing for a fixed $T_c$ (it emits a skyline only for `skyline`, `opt`, and `const`). The divergence is intentional and does not target v0 parity.

## Rationale

- The skyline optimizer produces per-segment $T_c$ values on explicit boundaries; per-segment rows carry that structure directly, without interpolating onto a date grid.
- Reporting $T_c$ alongside $N_e$ exposes the inferred time scale, the quantity the optimizer and the confidence band ([kb/decisions/coalescent-skyline-hessian-confidence-bands.md](coalescent-skyline-hessian-confidence-bands.md)) actually operate on.
- One rich document with a flat projection keeps JSON and delimited output consistent by construction.

## Limits

The schema is not interchangeable with v0's `skyline.tsv`. Consumers expecting v0's `date`/`N_e`/`lower`/`upper` columns must be updated to the per-segment layout.
