# Skyline plot output missing and grid-points default diverges from v0

Let $T_c$ denote the coalescent population-size time scale. The inferred coalescent time scale (skyline, optimized-constant, or fixed $T_c$) is serialized from the `timetree` command to `timetree.coalescent.tsv` (in the default `--output-all` set) and, on request, `timetree.coalescent.csv` and `timetree.coalescent.json`. Each carries the per-segment $T_c$, effective population size $N_e$, and the confidence band. Two gaps to v0 remain.

## Remaining gaps

- No plot. v0 `--coalescent=skyline` also renders `skyline.pdf`. v1 produces no plot; skyline plotting is part of the unimplemented plotting surface (see [N-timetree-plot-unimplemented.md](N-timetree-plot-unimplemented.md)).
- Grid-points default. v0 defaults to 20 skyline grid points; v1 defaults to 10 (`--skyline-n-points` in [packages/treetime/src/commands/timetree/args.rs](../../packages/treetime/src/commands/timetree/args.rs)). Whether v1 should adopt the v0 default is undecided.

## Notes

The serialized format is v1's own rich (JSON) and flat (TSV/CSV) layout with `segment.start` / `segment.end`, `T_c`, and `N_e` fields, not v0's `#date`, `N_e`, `lower`, `upper` columns. The optimizer and confidence band also diverge from v0, so skyline values differ (see [kb/decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md) and [kb/decisions/coalescent-skyline-hessian-confidence-bands.md](../decisions/coalescent-skyline-hessian-confidence-bands.md)).
