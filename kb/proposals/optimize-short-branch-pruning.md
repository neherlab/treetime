# Short branch pruning for optimize command

After the branch-length optimization loop converges, prune internal branches too short to be resolved by the available sequence data. Matches v0's marginal-mode behavior. Listed as `[ ]` in [kb/features/optimize.md](../features/optimize.md) (line 73).

## Criterion

v0 prunes an internal edge when both conditions hold:

- $b < 0.1m$, where $L$ is the sequence length and $m=1/L$ is the one-mutation resolution
- the zero-time sequence-transition likelihood for the reconstructed parent and child sequences is greater than $0.1$

The first condition identifies branches shorter than 10% of a single substitution. The second evaluates `def GTR.prob_t()` at branch length zero with the reconstructed sequences and their pattern multiplicities, matching `def TreeAnc.prune_short_branches()`. [`packages/legacy/treetime/treetime/treeanc.py#L1475-L1495`](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1495)

`fn is_zero_branch_optimal()` evaluates a derivative sign and is not equivalent to the v0 likelihood threshold. The pruning implementation must compute the v0 predicate directly.

## Pipeline insertion

After `run_optimize_loop` completes and before the final `update_marginal` pass:

1. Compute the one-mutation resolution from the sequence length.
2. For each eligible internal edge, evaluate both v0 pruning conditions.
3. Collapse marked edges (merge child into parent, creating polytomies)
4. Reconcile partition topology (`reconcile_topology`, same pattern as `timetree/refinement.rs`)
5. Run final `update_marginal` on the pruned tree

v0 marginal mode prunes after the loop. v0 joint mode (removed in v1) pruned inside the loop.

## Relation to existing zero-branch handling

The optimization loop already calls `find_zero_optimal_internal_edges` per iteration, which identifies edges where the optimizer drives the branch length to exactly zero. These are collapsed during the loop as part of topology cleanup. Post-loop pruning is complementary: it catches edges that converged to a small but nonzero length that the alignment cannot distinguish from zero.

## Locations

- Pipeline insertion: `packages/treetime/src/optimize/pipeline.rs`
- Probability calculation: use the same reconstructed states and pattern multiplicities as v0.
- Topology collapse: `fn collapse_edge()` [`packages/treetime/src/optimize/topology/collapse.rs#L34-L82`](../../packages/treetime/src/optimize/topology/collapse.rs#L34-L82)
- Partition reconciliation: `PartitionRerootOps::reconcile_topology` in `partition/traits.rs`

## Validation

- flu/h3n2/20: count pruned branches, compare with v0 marginal-mode output
- Synthetic tree with intentionally short branches: verify threshold matches `0.1 * one_mutation`
- All branches above threshold: no pruning, output identical to current
- The root has no incoming edge and is excluded; eligible internal children of the root use the ordinary predicate.

## Related

- [kb/proposals/optimize-pipeline-timetree-parity.md](optimize-pipeline-timetree-parity.md) -- parent proposal
- [kb/features/optimize.md](../features/optimize.md) -- `[ ] Short branch pruning after optimization` and `[ ] MIN_BRANCH_LENGTH floor for GTR calculations`
- [kb/tickets/optimize-add-short-branch-pruning.md](../tickets/optimize-add-short-branch-pruning.md) -- implementation ticket
