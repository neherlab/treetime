# Timetree plot outputs are unimplemented

`timetree` renders no plot, and it has no plot flags (v0: `--plot-rtt`, `--plot-tree`, `--tip-labels`, `--no-tip-labels`) until the plots exist.

A public enhancement request asks TreeTime to emit root-to-tip plots at several analysis points, so that date parsing, clock filtering, and failures can be diagnosed [[issue](https://github.com/neherlab/treetime/issues/228)].

## Missing plots

- **Root-to-tip regression** (`--plot-rtt`): the `clock` command already renders this chart as SVG, PNG, or PDF by the file extension of `--output-rtt-chart` (`RttChart` [packages/app-output/src/output_plan.rs#L715](../../packages/app-output/src/output_plan.rs#L715)), so `timetree` can reuse it
- **Time-scaled tree** (`--plot-tree`)
- **Skyline**: v0 `--coalescent=skyline` renders `skyline.pdf`. v1 writes the coalescent time scale $T_c$ (skyline or optimized constant) only as tables: `timetree.coalescent.tsv` in the default `--output-all` set of a run with a coalescent prior (a `.csv` file name writes commas), and on request `timetree.coalescent.json` [packages/app-output/src/output_plan.rs#L713-L714](../../packages/app-output/src/output_plan.rs#L713-L714). Each carries the per-segment $T_c$, the effective population size $N_e$, and the confidence band

The table layout is v1's own per-segment layout, not v0's `date` / `N_e` / `lower` / `upper` columns; [kb/decisions/coalescent-output-schema.md](../decisions/coalescent-output-schema.md) records the schema. The optimizer and the confidence band also differ from v0, so skyline values differ ([kb/decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md), [kb/decisions/coalescent-skyline-hessian-confidence-bands.md](../decisions/coalescent-skyline-hessian-confidence-bands.md)).

## Fix

Implement the three plots as output selections of `timetree`, and add their flags with them. Each plot draws through `packages/treetime-charts`, as the `clock` chart does: plotters draws the plot once as SVG with the embedded IBM Plex Sans font (`register_chart_fonts()`, `CHART_FONT_FAMILY`), with the colors of `palette`, the axis ticks of `ticks` and the text measurement of `text`, and `encode_chart()` in `packages/app-commands/src/rtt_chart/write.rs` turns that SVG into the format that the file extension of the output names (`ChartFormat` in `packages/app-output/src/output_format.rs`: SVG, PNG or PDF), so every plot gets the three formats.
