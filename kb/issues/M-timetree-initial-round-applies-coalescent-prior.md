# Timetree applies the coalescent prior before the refinement loop

When a coalescent model is set, v1's init phase infers the node times once under the mean-field prior seeded from the Gaussian node marginals, without a prior-free pass (`fn initialize` in [packages/treetime/src/timetree/init.rs](../../packages/treetime/src/timetree/init.rs)). v0's initial time tree runs without any coalescent prior, because the merger model is only set inside the loop ([packages/legacy/treetime/treetime/treetime.py#L270](../../packages/legacy/treetime/treetime/treetime.py#L270), [treetime.py#L1056](../../packages/legacy/treetime/treetime/treetime.py#L1056)).

[kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md) covers how lineage counts feed the prior, not when the prior is first applied.

## Impact

With any `--coalescent` prior (a fixed time scale, `opt` or `skyline`), the dates entering the first refinement round already include the prior, unlike v0.

The init phase is decided (phase 1 of
[timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md));
the divergence from v0 needs its decision file (phase 5), which closes this issue.
