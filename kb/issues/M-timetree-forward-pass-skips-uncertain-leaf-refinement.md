# Timetree forward pass skips uncertain leaf posterior refinement

The marginal forward pass returns early for every leaf node before combining the parent message with the leaf's time constraint. V0 skips only exact delta-constrained leaves; uncertain leaves (Range or Function distributions representing imprecise dates) receive a marginal posterior that incorporates information from the rest of the tree.

## Mechanism

In [`packages/treetime/src/timetree/inference/forward_pass.rs#L73`](../../packages/treetime/src/timetree/inference/forward_pass.rs#L73), `refine_distribution_from_parent` returns immediately for all leaves:

```rust
if graph.is_leaf(slot.key) {
    return Ok(());
}
```

This means no leaf -- regardless of whether its date constraint is an exact point or an uncertain range -- gets its distribution refined from the parent's cavity message. The leaf retains only its input constraint.

## V0 behavior

V0 `ClockTree._ml_t_marginal()` calls `_get_subtree_time_response()` for the forward pass. For each child, it computes the parent cavity distribution (parent divided by child's backward message), convolves with the branch length, and updates the child's time distribution. Leaves with exact (`Delta`) constraints retain their constraint; uncertain leaves receive the combined posterior.

## Impact

Uncertain leaves (imprecise dates, date ranges, missing dates with priors) do not benefit from tree-wide information. Their inferred dates are constrained only by the input, not by the marginal posterior that v0 computes. This is an unapproved parity defect affecting date inference quality for any dataset with uncertain sample dates.

## Related issues

- [M-inference-forward-backward-asymmetry.md](M-inference-forward-backward-asymmetry.md): other forward/backward asymmetries (normalization, graph mutation, BL clamping -- distinct mechanisms)
- [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md): topology violations from independent peak selection (orthogonal to leaf refinement)
