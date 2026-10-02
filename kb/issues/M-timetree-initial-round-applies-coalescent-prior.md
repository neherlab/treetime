# Timetree applies the coalescent prior before the refinement loop

When a coalescent model is set, v1 runs the initial time inference twice: once without the prior, then again with it (`fn run_initial_round` in [packages/treetime/src/timetree/round.rs](../../packages/treetime/src/timetree/round.rs), the `prior_wanted` block after the first `run_timetree`). v0's initial time tree runs without any coalescent prior, because the merger model is only set inside the loop ([packages/legacy/treetime/treetime/treetime.py#L270](../../packages/legacy/treetime/treetime/treetime.py#L270), [treetime.py#L1056](../../packages/legacy/treetime/treetime/treetime.py#L1056)).

[kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md) covers how lineage counts feed the prior, not when the prior is first applied.

## Impact

With `--coalescent`, `--coalescent-opt` or `--coalescent-skyline`, the dates entering the first refinement round already include the prior, unlike v0.
