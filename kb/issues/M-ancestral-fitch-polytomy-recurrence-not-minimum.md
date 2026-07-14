# Fitch recurrence is not minimum parsimony on multifurcations

The sparse Fitch backward pass intersects all child state sets and falls back to their union when the common intersection is empty [packages/treetime/src/ancestral/fitch_sub.rs#L38-L72](../../packages/treetime/src/ancestral/fitch_sub.rs#L38-L72). That recurrence is exact for bifurcations, but it can retain non-minimum states when a node has three or more children.

For child states $A$, $A$, and $C$, the intersection is empty and the union is $\{A,C\}$. Choosing parent state $A$ costs one change, whereas choosing $C$ costs two. The forward pass may nevertheless select $C$ from an ancestor, producing a non-minimum reconstruction and mutation placement. Sparse marginal initialization inherits the selected Fitch state.

V0 uses the same intersection-or-union recurrence [packages/legacy/treetime/treetime/treeanc.py#L638-L655](../../packages/legacy/treetime/treetime/treeanc.py#L638-L655). Correcting v1 therefore improves scientific correctness while intentionally diverging from the reference output, which requires explicit approval.

## Potential solutions

- O1. Use an exact finite-state Sankoff recurrence at nodes with three or more children and retain every minimum-cost parent state; keep the binary Fitch fast path.
- O2. Resolve multifurcations into bifurcations before Fitch reconstruction. This couples ancestral output to a topology transformation and can introduce arbitrary resolution choices.
- O3. Preserve v0 parity. This retains known non-minimum ancestral states and mutation placements.

## Recommendation

Approve O1 as an intentional v0 divergence. Validate it with the explicit $A,A,C$ counterexample, exhaustive finite-state scores on generated multifurcations, and a golden comparison that records the approved difference. No implementation ticket is ready until that scientific-output divergence is approved.

## Related issues

- [N-ancestral-fitch-site-classification-parallel-regression.md](N-ancestral-fitch-site-classification-parallel-regression.md)
- [N-test-coverage-gaps.md](N-test-coverage-gaps.md)
