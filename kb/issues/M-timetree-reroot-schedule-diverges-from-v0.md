# Timetree reroot and clock-filter schedule differs from v0

The order of rerooting, clock filtering, branch-length optimization and time inference before the refinement loop differs between v1 and v0. In v1, `fn run()` fits the first clock model in `fn estimate_initial_clock()` [packages/treetime/src/timetree/pipeline.rs#L279-L315](../../packages/treetime/src/timetree/pipeline.rs#L279-L315), then `fn run_pre_loop()` runs the other pre-loop steps [packages/treetime/src/timetree/pre_loop.rs#L31-L73](../../packages/treetime/src/timetree/pre_loop.rs#L31-L73). In v0, `TreeTime._run` (`packages/legacy/treetime/treetime/treetime.py`) runs all steps.

## Side by side

| #   | v0                                                                                                                                   | v1                                                              |
| --- | ------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------- |
| 1   | -                                                                                                                                    | clock estimate and reroot on the raw input branch lengths       |
| 2   | ML branch-length optimization, one iteration ([treetime.py#L243](../../packages/legacy/treetime/treetime/treetime.py#L243))          | ML branch-length pre-step                                       |
| 3   | reroot, no covariation ([treetime.py#L486-L489](../../packages/legacy/treetime/treetime/treetime.py#L486-L489))                      | reroot, no covariation                                          |
| 4   | clock filter                                                                                                                         | clock filter                                                    |
| 5   | reroot with covariation ([treetime.py#L512-L513](../../packages/legacy/treetime/treetime/treetime.py#L512-L513))                     | - (approved, see below)                                         |
| 6   | ML branch-length optimization ([treetime.py#L266](../../packages/legacy/treetime/treetime/treetime.py#L266))                         | ML branch-length post-step                                      |
| 7   | -                                                                                                                                    | reroot, outliers excluded (v1 has no covariation)              |
| 8   | time tree, no coalescent prior ([treetime.py#L270](../../packages/legacy/treetime/treetime/treetime.py#L270))                        | time tree, then again with the coalescent prior when one is set |
| 9   | reroot, ancestral reconstruction, time tree ([treetime.py#L278-L285](../../packages/legacy/treetime/treetime/treetime.py#L278-L285)) | -                                                               |

Row 5 is an approved difference: [kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md) removes covariation from the clock regression. No decision records the differences in the other rows. Row 8 is tracked separately in [M-timetree-initial-round-applies-coalescent-prior.md](M-timetree-initial-round-applies-coalescent-prior.md).

The table describes the date-based methods `least-squares` and `oldest`. The date-free methods (`min-dev`, `--outgroup`, `midpoint`) root once, in row 1, and the reroots of rows 3 and 7 are clock refits on that root for them ([kb/decisions/reroot-date-free-methods-root-once.md](../decisions/reroot-date-free-methods-root-once.md)).

## Impact

The root position and the branch lengths entering the first time inference can differ from v0, which changes every downstream estimate. With `--keep-root`, v1's first clock model is fitted on the raw input branch lengths and is reused by the clock filter and the first time inference.

> [!IMPORTANT]
> **Decision required.** Rows 1, 7 and 9 differ from v0 without an approved decision: v1 fits a clock and reroots on the raw input branch lengths before any optimization (row 1), reroots again after the post-step with the outliers excluded (row 7), and omits the v0 reroot, ancestral reconstruction and time tree that end the v0 pre-loop (row 9). The options for each row: restore the v0 step, or record the v1 step as an intentional change in `kb/decisions/`. Evidence: the v0 lines cited in the table, and the v1 order in `fn estimate_initial_clock()` and `fn run_pre_loop()`.
