# Skyline plot output missing

Let $T_c$ denote the coalescent population-size time scale. The inferred coalescent time scale (skyline or optimized-constant $T_c$) is serialized from the `timetree` command to `timetree.coalescent.tsv` (in the default `--output-all` set) and, on request, `timetree.coalescent.csv` and `timetree.coalescent.json`. Each carries the per-segment $T_c$, effective population size $N_e$, and the confidence band. One gap to v0 remains.

## Remaining gap

- No plot. v0 `--coalescent=skyline` also renders `skyline.pdf`. v1 produces no plot; skyline plotting is part of the unimplemented plotting surface (see [N-timetree-plot-unimplemented.md](N-timetree-plot-unimplemented.md)).

## Notes

The serialized format is v1's own per-segment layout, not v0's `date` / `N_e` / `lower` / `upper` columns; the schema is decided in [kb/decisions/coalescent-output-schema.md](../decisions/coalescent-output-schema.md). The `--skyline-n-points` default is 20, matching v0. The optimizer and confidence band also diverge from v0, so skyline values differ (see [kb/decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md) and [kb/decisions/coalescent-skyline-hessian-confidence-bands.md](../decisions/coalescent-skyline-hessian-confidence-bands.md)).
