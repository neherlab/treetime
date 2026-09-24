# Error suppression via unwrap_or_default and silent fallbacks

## Summary

Multiple locations suppress errors by substituting default values for semantically important data, or by silently proceeding after failures that should be reported.

## Instances

### unwrap_or_default() on semantically important data (8 instances)

- `clock/clock_filter.rs:41:` branch length defaults to 0.0
- `clock/reroot.rs:178:` branch length defaults to 0.0
- `clock/rtt.rs:36:` branch length defaults to 0.0
- `seq/div.rs:26:` parent divergence defaults to 0.0
- `seq/div.rs:28:` branch length defaults to 0.0
- `partition/marginal_discrete.rs:60:` node name defaults to empty string
- `optimize/topology/merge_shared_mutations.rs:211-212:` indels default to empty
- `timetree/optimization/relaxed_clock.rs:87:` coefficients default to zero

A default of 0.0 for branch length can produce division-by-zero downstream or silently exclude the branch from optimization. An empty node name makes the node invisible to output serialization.

### debug_assert_eq! stripped in release (compose_substitutions)

`packages/treetime/src/seq/mutation.rs:100:`

`debug_assert_eq!(ps.qry(), cs.reff(), ...)` is stripped in release builds. A broken substitution chain (where the query state of the parent substitution does not match the reference state of the child substitution) silently produces incorrect mutation annotations.

### branch_length().unwrap_or(one_mutation) silent fallback

`packages/treetime/src/timetree/inference/runner.rs:106:`

Edges with no branch length get `one_mutation` (= 1.0 / total_sites) as a fallback when building Poisson branch-length distributions. A missing branch length could indicate a tree-loading error or an uninitialized edge, but the fallback silently assigns a plausible value.

### infer_gtr_impl silently proceeds after non-convergence

`packages/treetime/src/gtr/infer_gtr/common.rs:157-158:`

Returns `Ok(...)` with a `warn!` log when GTR inference does not converge. No convergence flag in the result struct. Callers cannot distinguish a converged model from a non-converged one without parsing log output.

### Composition::adjust_count saturating_add_signed silently clamps

`packages/treetime/src/seq/composition.rs:76-79:`

Uses `saturating_add_signed` which silently clamps to 0 on underflow. A negative composition count indicates a data integrity problem that should be reported, not masked.

### Empty sequence for a missing node

- `partition/marginal/dense/partition.rs:161:` `extract_ancestral_sequence()` returns an empty sequence when the node has no marginal state
- `partition/marginal/sparse/partition.rs:57:` `node_sequence()` returns an empty sequence when the node has no marginal state

A missing node state is a broken invariant. An empty sequence is written to the output as if the node had no sites.

### Mass-sizing errors read as "not sizable"

`packages/treetime-distribution/src/distribution_ops/mass_domain.rs:40:`

`peak_normalized_if_mass_sizable()` discards a `total_mass()` error through `is_ok_and`, so a distribution with an undeclared tail or fewer than two grid points is treated the same as one without finite total mass. Callers then take the fallback window intended for the second case.

## Impact

- Silent data corruption in release builds from unchecked substitution chains
- Branch length 0.0 defaults cause division-by-zero or exclusion from optimization
- Non-converged GTR models used without caller awareness

## Related tickets

- [kb/tickets/safety-missing-convergence-failure-reporting.md](../tickets/safety-missing-convergence-failure-reporting.md)
- [kb/tickets/safety-unwrap-or-default-and-silent-fallbacks-on-important-data.md](../tickets/safety-unwrap-or-default-and-silent-fallbacks-on-important-data.md)
