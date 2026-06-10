# Short branch pruning for optimize command

After the branch-length optimization loop converges, prune internal branches too short to be resolved by the available sequence data. Matches v0's marginal-mode behavior. Listed as `[ ]` in [kb/features/optimize.md](../features/optimize.md) (line 73).

## Criterion

v0 prunes an internal edge when both conditions hold:

- `bl < 0.1 * one_mutation`, where `one_mutation = 1.0 / sequence_length`
- `P(zero) > 0.1`, where P(zero) is the probability that the true branch length is zero

The first condition identifies branches shorter than 10% of a single substitution -- too small for the alignment to resolve. The second condition checks whether the optimizer's likelihood surface supports zero length at that edge.

v1 already computes P(zero) via `is_zero_branch_optimal` in `optimize/zero_boundary.rs`, which evaluates the combined derivative at zero and checks whether the likelihood is concave there. The pruning criterion reuses this existing machinery.

## Pipeline insertion

After `run_optimize_loop` completes and before the final `update_marginal` pass:

1. Compute `one_mutation = 1.0 / sequence_length`
2. For each internal edge: check `bl < 0.1 * one_mutation` AND P(zero) > 0.1
3. Collapse marked edges (merge child into parent, creating polytomies)
4. Reconcile partition topology (`reconcile_topology`, same pattern as `timetree/refinement.rs`)
5. Run final `update_marginal` on the pruned tree

v0 marginal mode prunes after the loop. v0 joint mode (removed in v1) pruned inside the loop.

## Relation to existing zero-branch handling

The optimization loop already calls `find_zero_optimal_internal_edges` per iteration, which identifies edges where the optimizer drives the branch length to exactly zero. These are collapsed during the loop as part of topology cleanup. Post-loop pruning is complementary: it catches edges that converged to a small but nonzero length that the alignment cannot distinguish from zero.

## Locations

- Pipeline insertion: `packages/treetime/src/optimize/pipeline.rs`
- Zero-branch logic: `packages/treetime/src/optimize/zero_boundary.rs`
- Topology collapse: reuse `remove_node_if_trivial` from `treetime-graph/src/reroot.rs` or the existing collapse path in `optimize/run_loop.rs`
- Partition reconciliation: `PartitionRerootOps::reconcile_topology` in `partition/traits.rs`

## Validation

- flu/h3n2/20: count pruned branches, compare with v0 marginal-mode output
- Synthetic tree with intentionally short branches: verify threshold matches `0.1 * one_mutation`
- All branches above threshold: no pruning, output identical to current
- Root children below threshold: root edge must not be pruned

## Related

- [kb/proposals/optimize-pipeline-timetree-parity.md](optimize-pipeline-timetree-parity.md) -- parent proposal
- [kb/features/optimize.md](../features/optimize.md) -- `[ ] Short branch pruning after optimization` and `[ ] MIN_BRANCH_LENGTH floor for GTR calculations`
- [kb/tickets/optimize-add-short-branch-pruning.md](../tickets/optimize-add-short-branch-pruning.md) -- implementation ticket
