# Coalescent output uses a per-segment schema

v1 serializes the inferred coalescent time scale from the `timetree` command with its own per-segment schema. This differs from v0, which writes one interpolated row per date.

**Type**: Output-format change from v0.

**Status**: Approved by the maintainer team.

**v0**: `MergerModel` writes a single `skyline.tsv` with columns `date`, `N_e`, `lower`, `upper`, one row per interpolation date across the tree's span ([`packages/legacy/treetime/treetime/CLI_io.py#L248`](../../packages/legacy/treetime/treetime/CLI_io.py#L248)). The coalescent time scale $T_c$ is not reported; only the derived effective population size $N_e$ and its band appear.

**v1**: The serializers in [`packages/treetime/src/commands/timetree/output/coalescent.rs`](../../packages/treetime/src/commands/timetree/output/coalescent.rs) emit two shapes from one document:

- Rich JSON (`timetree.coalescent.json`): an `inputs` block (mode, generations per year, confidence width) and `outputs.segments`, each segment carrying a numeric-date `segment` interval and nested `T_c` and `N_e` objects with `value`/`lower`/`upper`.
- Flat delimited TSV/CSV (`timetree.coalescent.tsv`, `timetree.coalescent.csv`): one row per segment with dotted columns `segment.start`, `segment.end`, `T_c.value`, `T_c.lower`, `T_c.upper`, `N_e.value`, `N_e.lower`, `N_e.upper`. The TSV is in the default `--output-all` set; CSV and JSON are opt-in.

## Decision

Report the coalescent as per-segment intervals rather than per-date points, and report both $T_c$ and the derived $N_e = T_c \cdot \text{gen\_per\_year}$ with their confidence bands. The rich JSON is the single source of truth; the flat rows are its projection, shifting the 0-based JSON index to 1-based.

A fixed, user-supplied $T_c$ is not inferred and writes no coalescent output at all, matching v0 (which emits a skyline only for `skyline`, `opt`, and `const`) and v1's own stderr $N_e$ report.

## Rationale

- The skyline optimizer produces per-segment $T_c$ values on explicit boundaries; per-segment rows carry that structure directly, without interpolating onto a date grid.
- Reporting $T_c$ alongside $N_e$ exposes the inferred time scale, the quantity the optimizer and the confidence band ([kb/decisions/coalescent-skyline-hessian-confidence-bands.md](coalescent-skyline-hessian-confidence-bands.md)) actually operate on.
- One rich document with a flat projection keeps JSON and delimited output consistent by construction.

## Limits

The schema is not interchangeable with v0's `skyline.tsv`. Consumers expecting v0's `date`/`N_e`/`lower`/`upper` columns must be updated to the per-segment layout.
