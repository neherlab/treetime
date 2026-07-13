# Coalescent multiplication ordering diverges from v0 without empirical validation

The v0 and v1 backward passes apply the coalescent contribution at different points. Because the implementations also reconstruct grids differently, the two paths do not merely reassociate exact arithmetic: they evaluate and interpolate on different grids.

## Reference behavior

In v0 [`ClockTree._ml_t_marginal()`](../../packages/legacy/treetime/treetime/clock_tree.py), the backward pass:

1. collects all child messages;
2. multiplies them into `product_of_child_messages`;
3. evaluates `merger_model.node_contribution()` on that final product grid; and
4. multiplies the evaluated contribution into the child product.

In v1 [`propagate_distributions_backward_single_node()`](../../packages/treetime/src/timetree/inference/backward_pass.rs), `result` is seeded with the plain coalescent `Formula`, then each child message is multiplied into it sequentially. The first `Formula x Function` operation fixes a grid before the remaining children are known, and later products interpolate again.

## Unverified decision

[`kb/decisions/coalescent-multiplication-ordering.md`](../decisions/coalescent-multiplication-ordering.md) classifies this as the same mathematical result within floating-point precision and calls the interpolation error negligible because the coalescent contribution is smooth. The entry contains no v0/v1 empirical comparison demonstrating the claimed bound across binary trees, polytomies, Tc values, skyline histories, or grids with non-aligned knots.

The project parity rule treats an unverified divergence as a bug. Peak-relative negative-log multiplication preserves the likelihood shape, but it does not establish that the two application orders construct equivalent grids.

## Decision required

No implementation ticket is ready until the existing decision is explicitly retained, revised, or removed. Project rules require that choice to receive explicit user consent in the implementation session. If exact parity is selected, the ticket must prescribe child-first combination and v0-backed numerical validation. If the v1 order is retained, the decision needs empirical error bounds and a separate stable underflow design.

## Related issues

- [M-distribution-product-grid-resolution-diverges-from-v0.md](M-distribution-product-grid-resolution-diverges-from-v0.md)
