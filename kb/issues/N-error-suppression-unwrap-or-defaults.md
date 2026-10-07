# Error suppression via unwrap_or_default and silent fallbacks

## Summary

Multiple locations suppress errors by substituting default values for semantically important data, or by silently proceeding after failures that should be reported.

## Instances

### Defaults on semantically important data

- `branch_lengths.rs:4-16:` `branch_lengths_or_zero()` and `branch_length_or_zero()` map a missing branch length to 0.0. Every pipeline calls them before the marginal passes and before the divergence sums
- `clock/clock_filter.rs:37:` branch length defaults to 0.0 through `branch_length_or_zero()`
- `clock/pipeline.rs:68:` branch length defaults to 0.0 through `branch_length_or_zero()` in the divergences that the RTT output reads. The clock regression panics on the same input instead: [M-clock-regression-panics-on-missing-branch-length.md](M-clock-regression-panics-on-missing-branch-length.md)
- `partition/marginal/discrete/partition.rs:42:` and `partition/marginal/discrete/input.rs:30:` node name defaults to empty string
- `optimize/topology/merge_shared_mutations.rs:100:` and `:212:` indels default to an empty slice (`map_or(empty_indels, ..)`) when the edge has no observation
- `timetree/optimization/relaxed_clock.rs:30:` branch length defaults to 0.0
- `timetree/optimization/relaxed_clock.rs:63:` coefficients default to zero

A default of 0.0 for branch length can produce division-by-zero downstream or silently exclude the branch from optimization. An empty node name makes the node invisible to output serialization.

### debug_assert_eq! stripped in release (compose_substitutions)

`packages/treetime/src/seq/mutation.rs:314:`

`debug_assert_eq!(ps.qry(), cs.reff(), ...)` in `compose_substitutions()` is stripped in release builds. A broken substitution chain (where the query state of the parent substitution does not match the reference state of the child substitution) silently produces incorrect mutation annotations.

Fix: replace the `debug_assert_eq!` with a check that runs in all builds and returns an error. `compose_substitutions()` (`packages/treetime/src/seq/mutation.rs:286:`) already returns `Result`, so the error can propagate without an API change.

### branch_length().unwrap_or(one_mutation) silent fallback

`packages/treetime/src/timetree/inference/runner.rs:183:`

Edges with no branch length get `one_mutation` (= 1.0 / total_sites) as a fallback when building Poisson branch-length distributions. A missing branch length could indicate a tree-loading error or an uninitialized edge, but the fallback silently assigns a plausible value.

### infer_gtr_impl silently proceeds after non-convergence

`packages/treetime/src/gtr/infer_gtr/common.rs:70-74:`

Returns `Ok(...)` with a warning log (`progress_warn!`) when GTR inference does not converge. `struct InferGtrResult` (`packages/treetime/src/gtr/infer_gtr/common.rs:106:`) has no convergence flag. Callers cannot distinguish a converged model from a non-converged one without parsing log output.

Fix: add a convergence flag to `InferGtrResult` so callers can detect and handle non-convergence programmatically.

### Composition::adjust_count saturating_add_signed silently clamps

`packages/treetime/src/seq/composition.rs:70:`

Uses `saturating_add_signed` which silently clamps to 0 on underflow. A negative composition count indicates a data integrity problem that should be reported, not masked.

### Empty sequence for a missing node

- `partition/marginal/dense/partition.rs:170:` dense `extract_ancestral_sequence()` returns an empty sequence when the node has no marginal state
- `partition/marginal/sparse/partition.rs:109:` sparse `extract_ancestral_sequence()` does the same

A missing node state is a broken invariant. An empty sequence is written to the output as if the node had no sites.

### Mass-sizing errors read as "not sizable"

`packages/treetime-distribution/src/distribution_ops/mass_domain.rs:29-31:`

`peak_normalized_if_mass_sizable()` discards a `mass_profile()` error through `let Ok(..) else { return None; }`, so a distribution with an undeclared tail or fewer than two grid points is treated the same as one without finite total mass. Callers then take the fallback window intended for the second case.

## Impact

- Silent data corruption in release builds from unchecked substitution chains
- Branch length 0.0 defaults cause division-by-zero or exclusion from optimization
- Non-converged GTR models used without caller awareness
