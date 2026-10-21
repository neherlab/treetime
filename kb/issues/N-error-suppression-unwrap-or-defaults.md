# Error suppression via defaults and silent fallbacks

Several locations substitute a default value for data that carries meaning, or continue silently where the data is inconsistent. None of them is known to change an output today, but each one hides the cause when an upstream defect occurs.

## Defaults on data that carries meaning

- [packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L97](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L97) and [packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L207](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L207): indels default to an empty slice (`map_or(empty_indels, ..)`) when the edge has no observation. An edge without an observation has no mutations, so this default is plausibly correct; a type that cannot be absent would make that explicit
- [packages/treetime/src/timetree/optimization/relaxed_clock.rs#L30](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs#L30): a missing branch length becomes 0.0. The denominator is `opt_len + one_mutation`, so no division by zero follows, but the edge silently enters the relaxed-clock fit as a zero-length branch
- [packages/treetime/src/timetree/optimization/relaxed_clock.rs#L63](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs#L63): missing coefficients default to zero

## `branch_length().unwrap_or(one_mutation)` fallback

[packages/treetime/src/timetree/inference/runner.rs#L277](../../packages/treetime/src/timetree/inference/runner.rs#L277): an edge without a branch length gets `one_mutation` (1 / total sites) when the Poisson branch-length distribution is built. A missing branch length can indicate a tree-loading error or an uninitialized edge, but the fallback silently assigns a plausible value.

## `Composition::adjust_count` clamps negative counts

[packages/treetime/src/seq/composition.rs#L61](../../packages/treetime/src/seq/composition.rs#L61) uses `saturating_add_signed`, which silently clamps a count to 0 on underflow. A negative composition count indicates a data integrity problem that should be reported, not masked.

## Fix

Replace each default with an explicit error that names the node or edge, or with a type that makes the missing case impossible, unless the default is provably correct (then state the invariant in the type).

## Related issues

- [N-seq-substitution-chain-checked-only-in-debug-builds.md](N-seq-substitution-chain-checked-only-in-debug-builds.md)
- [N-gtr-inference-non-convergence-not-reported-to-callers.md](N-gtr-inference-non-convergence-not-reported-to-callers.md)
- [N-partition-missing-node-state-yields-empty-sequence.md](N-partition-missing-node-state-yields-empty-sequence.md)
