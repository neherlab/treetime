# Timetree has no joint time inference, so default dates differ from v0

TreeTime v0 infers node times jointly by default. v1 runs only marginal time inference and has no joint mode, so `treetime timetree` with default flags reports different dates than v0. No decision records this divergence.

## v0 behavior

- The CLI default is `--time-marginal false` ([packages/legacy/treetime/treetime/argument_parser.py#L258-L265](../../packages/legacy/treetime/treetime/argument_parser.py#L258-L265)): "For 'false' or 'never', TreeTime uses the jointly most likely values for the divergence times"
- `never`, `only-final` and `confidence-only` run joint time inference in every round ([packages/legacy/treetime/treetime/treetime.py#L218-L222](../../packages/legacy/treetime/treetime/treetime.py#L218-L222)); `make_time_tree` calls `_ml_t_joint` unless marginal is requested ([packages/legacy/treetime/treetime/clock_tree.py#L420-L423](../../packages/legacy/treetime/treetime/clock_tree.py#L420-L423))
- `only-final` and `confidence-only` add one final marginal round after the loop ([packages/legacy/treetime/treetime/treetime.py#L390-L396](../../packages/legacy/treetime/treetime/treetime.py#L390-L396)). That round re-estimates the clock and overwrites the dates with the marginal peaks

## v1 behavior

- `run_timetree` always runs sum-product message passing with cavity division ([packages/treetime/src/timetree/inference/forward_pass.rs](../../packages/treetime/src/timetree/inference/forward_pass.rs)), and the node time is the marginal peak clamped to the parent time
- `TimeMarginalMode` only controls an extra pass and whether confidence intervals are extracted ([packages/treetime/src/timetree/pipeline.rs](../../packages/treetime/src/timetree/pipeline.rs), the `TimeMarginalMode::OnlyFinal` block and `extract_confidence_intervals`). `never` and `always` therefore give identical dates
- `only-final` runs an extra `run_timetree`, commits clock branch lengths undamped, and runs `marginal_update_timetree`; v0's final round has no ancestral update

## Impact

- Default timetree dates differ from v0 on every dataset
- `--time-marginal` has no effect on dates except through the extra `only-final` pass, whose semantics differ from v0

## Related KB entries

- [kb/decisions/ancestral-joint-reconstruction-removed.md](../decisions/ancestral-joint-reconstruction-removed.md) covers ancestral sequences only, not node times
- [kb/features/timetree.md](../features/timetree.md) marks `never` as done "(joint most-likely times)", which the code does not implement

## Open question

Either implement joint time inference (v0 `_ml_t_joint`, [packages/legacy/treetime/treetime/clock_tree.py](../../packages/legacy/treetime/treetime/clock_tree.py)) and give `only-final` the v0 meaning, or approve marginal-only time inference in `kb/decisions/` and define `only-final` for a marginal-only engine.
