# Coalescent multiplication ordering

This document describes differences in how v0 and v1 combine child time distributions with coalescent prior contributions during the backward pass of timetree inference.

**Type**: Implementation difference (same mathematical result within floating-point precision).

**v0 location**: `_ml_t_marginal()` in `packages/legacy/treetime/treetime/clock_tree.py:684-785`.

**v1 location**: `propagate_distributions_backward_single_node()` in `packages/treetime/src/commands/timetree/inference/backward_pass.rs:34-86`.

**Affected datasets**: All datasets with `coalescent_tc` enabled.

## Background: coalescent priors in timetree inference

Timetree inference estimates divergence times (internal node ages) from a phylogenetic tree with known topology and branch lengths. The algorithm propagates time distributions from leaves toward the root (backward pass), then from root toward leaves (forward pass). At each internal node, time constraints from all children must be combined. TreeTime implements this as maximum-likelihood phylodynamic analysis [Sagulenko et al. 2018].

The Kingman coalescent [Kingman 1982] provides a prior probability over tree node ages based on population genetics. Under the coalescent model, lineages merge backward in time according to a stochastic process governed by the effective population size. The model assigns higher probability to node times that are consistent with expected coalescence waiting times.

For k lineages at time t, the total coalescence rate (rate at which any pair merges) is [Kingman 1982, Wakeley 2009]:

```
λ(t) = k(t) · (k(t) - 1) / (2 · Tc(t))
```

where Tc(t) is the effective population size (coalescence time scale). The integral merger rate I(t) accumulates the probability of no coalescence up to time t:

```
I(t) = ∫₀ᵗ κ(t') dt'
```

where κ(t) = (k(t) - 1) / (2 · Tc(t)) is the per-lineage merger rate (rate at which one specific lineage merges with any other).

At an internal node with m children, the coalescent contribution to the node's time distribution is:

- **Probability space**: λ(t)^(m-1) · exp(-I(t) · (m-1))
- **Neg-log space**: (m-1) · (I(t) - log(λ(t)))

This contribution is multiplied with the product of child messages to produce the node's time distribution.

## Multiplication ordering

Both implementations compute the same mathematical product. The difference is the order of operations and the resulting grid points.

**v0** multiplies all child messages first, then multiplies the coalescent contribution into the product:

1. Collect child messages: `msgs_to_multiply = [child.marginal_pos_Lx for child in node.clades]` (line 714-716)
2. Multiply children: `product_of_child_messages = Distribution.multiply(msgs_to_multiply)` (line 726)
3. Evaluate coalescent on product's grid: `time_points = node.product_of_child_messages.x` (line 736)
4. Compute contribution: `merger_contribution = self.merger_model.node_contribution(node, time_points)` (line 743)
5. Final multiply: `subtree_distribution = Distribution.multiply([merger_contribution, product_of_child_messages])` (line 744-746)

**v1** starts with the coalescent contribution and multiplies child messages sequentially:

1. Initialize with coalescent: `result = coalescent_contribs.get(&node.key).to_plain()` (line 48-51)
2. Loop over children, multiplying each message into result: `result = distribution_multiplication(&current, &parent_message_arc)` (line 72-73)

## Grid handling differences

The implementations differ in how they construct the output grid for multiplication.

**v0** (`Distribution.multiply` in `packages/legacy/treetime/treetime/distribution.py:82-149`) uses the union of all input grids:

```python
x_vals = np.unique(np.concatenate([k.x for k in dists]))
```

All knot points from all input distributions are preserved. The output contains every unique x-value from every input.

**v1** (`multiply_function_function` in `packages/treetime-distribution/src/distribution_ops/multiply.rs:138-166`) uses a uniform grid:

```rust
let n_points = a.t().len().max(b.t().len());
let t = overlap_min + (overlap_max - overlap_min) * (i as f64 / (n_points - 1) as f64);
```

The output grid has max(len(a), len(b)) evenly-spaced points over the intersection of supports. Input distributions are interpolated to these grid points.

## Coalescent representation

**v0** uses eager evaluation. `Coalescent.node_contribution(node, time_points)` in `packages/legacy/treetime/treetime/merger_models.py:239-250` evaluates the coalescent formula on the passed grid and returns a `NodeInterpolator` (discretized array):

```python
y = (self.integral_merger_rate(t) - np.log(self.total_merger_rate(t))) * multiplicity
return NodeInterpolator(t, y, is_log=True)
```

The grid comes from `product_of_child_messages.x`, so the coalescent is evaluated at all the knot points that resulted from combining child messages.

**v1** uses lazy evaluation. `compute_internal_contribution_single` in `packages/treetime/src/commands/timetree/coalescent/contributions.rs:93-146` returns a `Distribution::Formula` - a closure that captures the merger rate interpolators and evaluates the formula on demand:

```rust
let eval_fn = move |t: f64| -> eyre::Result<f64> {
    let i_t = integral_merger_rate.eval(t);
    let log_lambda_t = lambda_t[0].ln();
    let neg_log_contrib = multiplicity * (i_t - log_lambda_t);
    Ok(neg_log_contrib)
};
Ok(Distribution::Formula(DistributionFormula::new(eval_fn, t_min, t_max)))
```

When this Formula meets a Function during multiplication (`multiply_formula_function` in `packages/treetime-distribution/src/distribution_ops/multiply.rs:193-221`), the Formula is evaluated at each of the Function's grid points.

## Numerical impact

Multiplication is commutative and associative. Both orderings compute the same mathematical product: `coalescent(t) · child₁(t) · child₂(t) · ... · childₙ(t)`.

The numerical path differs in interpolation:

- **v0**: One interpolation step at the final multiply. The coalescent is evaluated directly on the union grid, preserving all knot points.
- **v1**: Each multiplication produces a new uniform grid. The coalescent is evaluated on the first child's grid. Subsequent children introduce linear interpolation at each step.

For binary nodes (2 children, the common case in phylogenetic trees), both orderings perform similar interpolation. For polytomies (3+ children), v1 accumulates one extra interpolation per child.

The coalescent contribution is a smooth function - the integral I(t) and log(λ(t)) are both derived from piecewise-linear merger rates. Smooth functions are insensitive to grid density, so interpolation error is negligible on the 200+ point grids used in practice.

## Why the difference exists

v1's backward pass processes children in a single loop (lines 53-77 in `backward_pass.rs`), interleaving edge message storage (`set_msg_to_parent`) with accumulation into `result`. Separating these into two passes (collect all messages, then reduce) would match v0's ordering but adds complexity with no measurable benefit.

The v1 grid strategy (uniform grid with max points) was chosen for simplicity and predictable memory usage. The v0 strategy (union of grids) can produce arbitrarily large output grids when many distributions are multiplied, though this is mitigated by the grid thinning logic at line 107-119 in `distribution.py`.

## References

- Kingman, J.F.C. (1982). "The coalescent." Stochastic Processes and their Applications 13(3): 235-248. doi:10.1016/0304-4149(82)90011-4
- Sagulenko, P., Puller, V., Neher, R.A. (2018). "TreeTime: Maximum-likelihood phylodynamic analysis." Virus Evolution 4(1): vex042. doi:10.1093/ve/vex042
- Wakeley, J. (2009). "Coalescent Theory: An Introduction." Roberts & Company Publishers. ISBN 978-0-9747077-5-4
- Birkner, M., Blath, J. (2014). "Coalescence 2.0: a multiple branching of recent developments." arXiv:1401.5248
