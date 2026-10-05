# Timetree infers node times by marginal inference only

## Deviation

v1 infers node times only by marginal inference: sum-product message passing over the tree, with the node time at the peak of its marginal posterior, clamped to the parent time. v0 infers node times jointly by default (the jointly most likely times) and uses marginal inference only on request. v1 has no joint time inference, and none is planned.

## v0 behavior

- The CLI default is `--time-marginal false` ([packages/legacy/treetime/treetime/argument_parser.py#L258-L265](../../packages/legacy/treetime/treetime/argument_parser.py#L258-L265)): "For 'false' or 'never', TreeTime uses the jointly most likely values for the divergence times"
- `never`, `only-final` and `confidence-only` run joint time inference in every round ([packages/legacy/treetime/treetime/treetime.py#L218-L222](../../packages/legacy/treetime/treetime/treetime.py#L218-L222)); `make_time_tree` calls `_ml_t_joint` unless marginal inference is requested ([packages/legacy/treetime/treetime/clock_tree.py#L420-L423](../../packages/legacy/treetime/treetime/clock_tree.py#L420-L423))
- `only-final` and `confidence-only` add one final marginal round after the loop ([packages/legacy/treetime/treetime/treetime.py#L390-L396](../../packages/legacy/treetime/treetime/treetime.py#L390-L396)). That round re-estimates the clock and overwrites the dates with the marginal peaks

## v1 behavior

- `run_timetree` runs the backward and forward sum-product passes in every round ([packages/treetime/src/timetree/inference/forward_pass.rs](../../packages/treetime/src/timetree/inference/forward_pass.rs), [packages/treetime/src/timetree/inference/backward_pass.rs](../../packages/treetime/src/timetree/inference/backward_pass.rs)), whatever `--time-marginal` says
- `--time-marginal` controls only the extra final round of `only-final` and whether confidence intervals are extracted. What each value should mean with marginal-only inference is open: [kb/issues/M-timetree-time-marginal-modes-undefined-for-marginal-only-inference.md](../issues/M-timetree-time-marginal-modes-undefined-for-marginal-only-inference.md)

## Rationale

v1 keeps one time-inference engine. The sum-product passes also produce the posterior distribution of every node time, which confidence-interval extraction reads, so the marginal engine is needed in every configuration. Joint inference would be a second engine (max-product passes with a traceback) that serves only the default dates.

## Consequences

- Default `timetree` dates differ from v0 on every dataset, because v0's default is joint inference. Comparisons with v0 use `--time-marginal always` on the v0 side
- Inputs exist where v0's default succeeds and marginal inference fails. On `data/mpox/clade-ii/20` with a fixed clock rate of 0.5 (about 5000 times a realistic mpox rate) and one iteration, v0's default finishes, and v0 with `--time-marginal always` fails with "Unexpected behavior detected in multiply function when determining peak of function with y-values '[]'". v1 fails on the same input
- [kb/decisions/ancestral-joint-reconstruction-removed.md](ancestral-joint-reconstruction-removed.md) covers joint ancestral sequence reconstruction; this entry covers node times
