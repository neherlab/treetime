# Log-likelihood values are indistinguishable from raw floats

Aggregate log-likelihood values use `f64` across partition, coalescent, optimization, mugration, and timetree convergence APIs. The compiler therefore accepts accidental mixing with probabilities, linear likelihoods, branch lengths, dates, and unrelated optimization scalars.

Representative boundaries include:

- `HasLogLh::get_log_lh()` and `graph_log_lh()` in [`packages/treetime/src/partition/traits.rs`](../../packages/treetime/src/partition/traits.rs)
- `DenseSeqDistribution::log_lh` in [`packages/treetime/src/partition/dense.rs`](../../packages/treetime/src/partition/dense.rs)
- `SparseSeqDistribution::log_lh` in [`packages/treetime/src/partition/sparse.rs`](../../packages/treetime/src/partition/sparse.rs)
- coalescent totals in [`packages/treetime/src/coalescent/edge_data.rs`](../../packages/treetime/src/coalescent/edge_data.rs) and [`packages/treetime/src/coalescent/total_lh.rs`](../../packages/treetime/src/coalescent/total_lh.rs)
- convergence components in [`packages/treetime/src/timetree/convergence/metrics.rs`](../../packages/treetime/src/timetree/convergence/metrics.rs)

Explicit `log_lh` naming reduces ambiguity for readers but cannot prevent a raw value from crossing the wrong API boundary. Such mistakes can preserve valid floating-point ranges and produce plausible scientific output, so ordinary runtime failures are not a reliable safeguard.

## Required resolution

Introduce the transparent `LogLh` scalar described in [kb/proposals/typed-log-likelihood-values.md](../proposals/typed-log-likelihood-values.md) and propagate it through aggregate scalar likelihood boundaries. Keep probability arrays and numerical kernels on native `f64`/`ndarray` representations.

The change must preserve numerical results, serialized scalar representation, and runtime performance. It must not add implicit mixed arithmetic, global finite/non-positive restrictions, or compatibility aliases.

## Related ticket

- [kb/tickets/likelihood-add-log-lh-scalar-type.md](../tickets/likelihood-add-log-lh-scalar-type.md)
