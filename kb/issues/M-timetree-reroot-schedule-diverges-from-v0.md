# Timetree reroot and clock-filter schedule differs from v0

The order of rerooting, clock filtering, branch-length optimization and time inference before the refinement loop differs between v1 (`packages/treetime/src/timetree/pipeline.rs`, fn `run`) and v0 (`packages/legacy/treetime/treetime/treetime.py`, `TreeTime._run`). No decision records the differences.

## Side by side

| #   | v0                                                                                                                                   | v1                                                              |
| --- | ------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------- |
| 1   | -                                                                                                                                    | clock estimate and reroot on the raw input branch lengths       |
| 2   | ML branch-length optimization, one iteration ([treetime.py#L243](../../packages/legacy/treetime/treetime/treetime.py#L243))          | ML branch-length pre-step                                       |
| 3   | reroot, no covariation ([treetime.py#L486-L489](../../packages/legacy/treetime/treetime/treetime.py#L486-L489))                      | reroot, no covariation                                          |
| 4   | clock filter                                                                                                                         | clock filter                                                    |
| 5   | reroot with covariation ([treetime.py#L512-L513](../../packages/legacy/treetime/treetime/treetime.py#L512-L513))                     | -                                                               |
| 6   | ML branch-length optimization ([treetime.py#L266](../../packages/legacy/treetime/treetime/treetime.py#L266))                         | ML branch-length post-step                                      |
| 7   | -                                                                                                                                    | reroot with covariation                                         |
| 8   | time tree, no coalescent prior ([treetime.py#L270](../../packages/legacy/treetime/treetime/treetime.py#L270))                        | time tree, then again with the coalescent prior when one is set |
| 9   | reroot, ancestral reconstruction, time tree ([treetime.py#L278-L285](../../packages/legacy/treetime/treetime/treetime.py#L278-L285)) | -                                                               |

Row 8 is tracked separately in [M-timetree-initial-round-applies-coalescent-prior.md](M-timetree-initial-round-applies-coalescent-prior.md).

## Impact

The root position and the branch lengths entering the first time inference can differ from v0, which changes every downstream estimate. With `--keep-root`, v1's first clock model is fitted on the raw input branch lengths and is reused by the clock filter and the first time inference.

## Related issues

- [M-clock-reroot-policy-boolean-selects-workflows.md](M-clock-reroot-policy-boolean-selects-workflows.md)
