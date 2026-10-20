# Timetree plot outputs are unimplemented

`timetree` renders no plot. The hidden flags `--plot-rtt` and `--plot-tree` are accepted by clap [packages/app-commands/src/commands/timetree/args.rs#L299-L311](../../packages/app-commands/src/commands/timetree/args.rs#L299-L311) and then rejected during argument conversion with "not yet implemented" through `make_error!()` [packages/app-commands/src/commands/timetree/args.rs#L58-L63](../../packages/app-commands/src/commands/timetree/args.rs#L58-L63).

A public enhancement request asks TreeTime to emit root-to-tip plots at several analysis points, so that date parsing, clock filtering, and failures can be diagnosed [[issue](https://github.com/neherlab/treetime/issues/228)].

## Missing plots

- **Root-to-tip regression** (`--plot-rtt`): the `clock` command already renders this chart as SVG and PNG (`ClockChartSvg`, `ClockChartPng` [packages/app-output/src/output_plan.rs#L846-L847](../../packages/app-output/src/output_plan.rs#L846-L847)), so `timetree` can reuse it
- **Time-scaled tree** (`--plot-tree`)
- **Skyline**: v0 `--coalescent=skyline` renders `skyline.pdf`. v1 writes the coalescent time scale $T_c$ (skyline or optimized constant) only as tables: `timetree.coalescent.tsv` in the default `--output-all` set, and on request `timetree.coalescent.csv` and `timetree.coalescent.json` [packages/app-output/src/output_plan.rs#L843-L845](../../packages/app-output/src/output_plan.rs#L843-L845). Each carries the per-segment $T_c$, the effective population size $N_e$, and the confidence band

The table layout is v1's own per-segment layout, not v0's `date` / `N_e` / `lower` / `upper` columns; [kb/decisions/coalescent-output-schema.md](../decisions/coalescent-output-schema.md) records the schema. The optimizer and the confidence band also differ from v0, so skyline values differ ([kb/decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md), [kb/decisions/coalescent-skyline-hessian-confidence-bands.md](../decisions/coalescent-skyline-hessian-confidence-bands.md)).

## Fix

Implement the three plots as output selections of `timetree` with the `plotters` owner that the `clock` chart uses, and remove the "not yet implemented" rejection, or remove the two hidden flags until the plots exist.
