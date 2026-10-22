# Polytomy mutation counts read the unconstrained branch length

Timetree edges carry two branch lengths: `branch_length`, the free ML or input estimate, and
`clock_branch_length`, the length the inferred node times imply. See
[timetree-clock-constrained-profile-propagation.md](../decisions/timetree-clock-constrained-profile-propagation.md).
Polytomy resolution wants the clock-constrained value but reads the free one.

## Mutation counts fall back to the free length

`fn edge_mutation_count()`
[packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L198-L222](../../packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L198-L222)
estimates a branch's substitution count as `round(branch_length * total_length)` when the
reconstructed substitution list cannot be read. The caller passes the free `branch_length`
([packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L176-L186](../../packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L176-L186)).
The estimate should use the clock-constrained length for consistency with the sweep it feeds,
which places those substitutions in time against a coalescent rate defined on the same time axis.

The fallback runs on every edge in input branch-length mode (`BranchModel::Input` has no
reconstruction), and in marginal mode on an edge for which a sparse partition stores no
substitution estimates.

## Not affected

The transition counts of GTR inference (`fn count_transitions_dense()` in
[packages/treetime/src/partition/marginal/shared/data.rs#L15](../../packages/treetime/src/partition/marginal/shared/data.rs#L15))
also read the ML length. That is unobservable: GTR inference runs only before the refinement loop,
so the two lengths are equal at every point where the counts are used. It becomes a real choice only
if GTR re-estimation moves into the loop.

## Impact

Polytomy resolution mis-weights branches in the fallback path, which input branch-length mode
takes on every edge.

## Validation

The change alters the sampled resolutions and therefore node times; it needs a before/after
comparison of input branch-length mode on a dataset with polytomies, such as
`data/mpox/clade-ii/1000`, besides a unit test.
