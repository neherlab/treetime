# Coalescent multiplication ordering

| Property    | Value                                                                                                                                                                                                                                                  |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Implementation difference (same results within floating-point precision)                                                                                                                                                                               |
| v1 Location | `propagate_distributions_backward_single_node()` (`#propagate_distributions_backward_single_node`) in [`packages/treetime/src/commands/timetree/inference/backward_pass.rs`](../../packages/treetime/src/commands/timetree/inference/backward_pass.rs) |
| v0 Location | `_ml_t_marginal()` (`#_ml_t_marginal`) in [`packages/legacy/treetime/treetime/clock_tree.py#L684-L785`](../../packages/legacy/treetime/treetime/clock_tree.py#L684-L785)                                                                               |
| Affects     | Grid points used for coalescent evaluation, interpolation path                                                                                                                                                                                         |
| Datasets    | All datasets with `coalescent_tc` enabled                                                                                                                                                                                                              |

## Multiplication ordering

v0 multiplies all child messages first, then multiplies the coalescent
contribution into the product:

```
product = child1 * child2 * ... * childN     (on union of all child grids)
result  = product * coalescent               (coalescent evaluated on product's grid)
```

v1 starts with the coalescent contribution and multiplies child messages
sequentially:

```
result  = coalescent * child1                (coalescent evaluated on child1's grid)
result  = result * child2                    (Function * Function)
...
result  = result * childN                    (Function * Function)
```

## Why the difference exists

v1's backward pass processes children in a single loop, interleaving edge message
storage (`set_msg_to_parent`) with accumulation into `result`. Separating these
into two passes (collect all messages, then reduce) would match v0's ordering but
adds complexity with no measurable benefit.

## Numerical impact

Multiplication is commutative and associative. Both orderings compute the same
mathematical product `coalescent(t) * child1(t) * ... * childN(t)`.

The numerical path differs in grid selection and interpolation:

- v0 evaluates the coalescent on the union of all children's grids (all knot
  points preserved). One interpolation step for the final multiply.
- v1 evaluates the coalescent on the first child's grid. Each subsequent child
  multiplication introduces one linear interpolation step.

For binary nodes (2 children, the common case in phylogenetic trees), both
orderings perform the same number of interpolation steps. For polytomies (3+
children), v1 accumulates one extra interpolation per child. The coalescent
contribution is a smooth function (integral and log of piecewise-linear merger
rates), making it insensitive to grid density. The accumulated interpolation
error on 200+ point grids is negligible.

## Coalescent representation

v0 uses eager evaluation: `Coalescent.node_contribution(node, time_points)`
evaluates the coalescent formula on a passed grid and returns a `NodeInterpolator`
(discretized array). The grid comes from `product_of_child_messages.x`.

v1 uses lazy evaluation: `compute_coalescent_contributions()` returns
`Distribution::Formula` closures that capture merger rate interpolators. The
Formula is evaluated when it meets a `Function` during multiplication
(`multiply_formula_function` evaluates the Formula at each of the Function's grid
points).

Both approaches produce the same values at the same grid points. The difference
is when evaluation happens (creation time vs multiplication time).
