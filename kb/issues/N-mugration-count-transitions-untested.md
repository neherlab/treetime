# Mugration transition counting lacks a discrete hand oracle

Mugration reuses the marginal passes and the transition counting of the sequence partitions. The shared code is alphabet-agnostic today, but no test computes a complete discrete result by hand, and no test runs the shared code with a state count other than 2 or 4. A later change that assumes four states (a hardcoded alphabet size, an ACGT or gap assumption) could silently produce wrong discrete-trait posteriors or wrong discrete GTR estimates.

This is a coverage and hardening gap, not a demonstrated defect.

## Current state

Transition counting drives the GTR estimation of mugration. `fn count_transitions_dense()` [packages/treetime/src/partition/marginal/shared/data.rs#L15-L57](../../packages/treetime/src/partition/marginal/shared/data.rs#L15-L57) accumulates:

- `nij`: the expected transition count matrix from `get_branch_mutation_matrix` over all edges
- `Ti`: the dwell times per state
- `root_state`: the argmax of the root posterior per row, skipped when the row is near uniform and `filter_uninformative_root` is set

The dense and discrete partitions reach it through one blanket implementation, `impl<P: IndexedPartition> MarginalPasses for P`, whose `fn count_transitions()` calls it [packages/treetime/src/partition/marginal/shared/pass.rs#L29-L85](../../packages/treetime/src/partition/marginal/shared/pass.rs#L29-L85).

The shared code takes the state count from the model and the data, not from a nucleotide constant:

- The count sizes `nij` and `Ti` from `gtr.pi().len()` [packages/treetime/src/partition/marginal/shared/data.rs#L24](../../packages/treetime/src/partition/marginal/shared/data.rs#L24)
- The normalization of a degenerate row uses `1.0 / n_cols` [packages/treetime/src/partition/marginal/shared/normalize.rs#L16](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L16)
- The root-uniform filter threshold is `1/n_states + 1e-10` ([kb/decisions/mugration-root-state-filtering.md](../decisions/mugration-root-state-filtering.md))
- The discrete partition builds its leaf profiles with `one_hot_profile(index, n_states)` and `missing_trait_profile(n_states)` [packages/treetime/src/partition/marginal/discrete/input.rs#L10-L18](../../packages/treetime/src/partition/marginal/discrete/input.rs#L10-L18)

## Coverage gap

- The shared contract tests cover parent/child orientation, edge accumulation, branch-length scaling, root composition, diagonal zeroing, and dense/sparse equality [packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs#L74-L119](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs#L74-L119) [packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs#L413-L461](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs#L413-L461)
- `packages/treetime/src/mugration/__tests__/test_discrete_marginal.rs` runs the marginal passes on a 2-state alphabet only and does not call `count_transitions()`
- The `test_gm_mugration_outputs` golden masters and the `test_run_mugration_*` tests in `packages/treetime/src/mugration/__tests__/test_run.rs` cover the path only indirectly

## Proposed test

Build a 3-leaf tree with known traits and branch lengths, run the marginal passes through the discrete partition, call `count_transitions()`, and compare `nij`, `Ti`, and `root_state` with values computed by hand from the GTR transition equations. Run the same fixture with 2, 3, and 5 states, and compare the node posteriors as well.

## Acceptance criteria

- Hand-computed fixtures with 2, 3, and 5 states verify the full transition matrix, the dwell-time vector, the root state, and the node posteriors
- The hand oracle distinguishes near-uniform root filtering; the existing tests remain the oracle for orientation, diagonal zeroing, branch-length weighting, and dense/sparse equality
- The oracle does not call the production accumulation helpers under test

## Related

- [N-discrete-missing-zero-states-inf.md](N-discrete-missing-zero-states-inf.md)
- [kb/algo/mugration.md](../algo/mugration.md)
