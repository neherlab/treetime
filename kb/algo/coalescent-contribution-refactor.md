# Coalescent Contribution Refactor

This document describes the problems with the current implementation of
coalescent prior contributions and specifies a cleaner architecture.

---

## Problems with the current implementation

### One closure per node, all identical except for multiplicity

`compute_node_contributions` allocates a `DistributionFormula` closure for
every node in the tree. Each closure captures cloned `Arc` references to
`integral_merger_rate`, `lineage_counts`, and `tc_dist` — the same three
objects for every node. The only per-node variation is `multiplicity = n_children - 1`,
a single scalar. This produces O(N) heap allocations and O(N) Arc reference
counts for data that could be shared as a single struct.

### TBP coordinate conversion on every evaluation

The underlying functions (`integral_merger_rate`, `lineage_counts`) are
expressed in time-before-present (TBP) coordinates, but the backward pass
operates in calendar time. Every closure evaluation performs a
`CalendarTime → Tbp` conversion (`t_tbp = present_time - t_calendar`). This
is a single subtraction that could be baked in once at construction time by
re-expressing both piecewise functions in calendar time before storing them.

### Per-evaluation heap allocations in `compute_merger_rates`

The internal-node closure calls:
```rust
compute_merger_rates(&Array1::from_vec(vec![k_t]), &Array1::from_vec(vec![tc_t]))
```

This allocates two `Vec`s and two `Array1`s to compute what are scalar
arithmetic operations. Multiplied across 1000 grid points per node and N
nodes, this produces millions of unnecessary small allocations.

### Formula applied at the wrong point in the backward pass

The current backward pass seeds `result` with the coalescent contribution
Formula *before* processing children:

```
result = Some(coalescent_formula)
for child in children:
    result = result * child_message      # Formula × Function on first child
```

`Formula × Function` must discretize the Formula onto a grid. At this point
only one child's message is known; the grid chosen (`b.len()` points over the
overlap) may not match what the fully-combined message will look like. The
correct moment to apply the coalescent contribution is *after* all child
messages have been accumulated, when the result's grid is fully determined.

### Leaf contributions are computed but never applied

`compute_node_contributions` builds a contribution for every node including
leaves. But the backward pass unconditionally sets `result = None` for leaf
nodes and never retrieves their contribution from the map, so the leaf
coalescent contribution is silently discarded.

This is incorrect. Each branch contributes `exp(I(t_child) - I(t_parent))` to
the coalescent likelihood. Summing over all branches, every node's `exp(I(t))`
term appears with a net coefficient equal to `n_children - 1` for internal
nodes and `-1` for leaves (from the parent branch denominator). A leaf with
known date still contributes `exp(-I(t_leaf))` that must cancel the
`exp(I(t_leaf))` numerator from the branch above it. Omitting leaf terms
breaks the cancellation and produces an incorrect likelihood.

---

## Proposed architecture

### Replace N closures with one evaluable struct

Define a `CoalescentModel` struct that holds the precomputed functions in
**calendar time** and exposes a scalar evaluation method:

```
struct CoalescentModel {
    // I(t) re-expressed in calendar time
    integral_merger_rate: PiecewiseLinearFn,
    // k(t) re-expressed in calendar time
    lineage_counts: PiecewiseConstantFn,
    // Tc(t) — constant or skyline, calendar time
    tc: Distribution,
    // domain: [earliest_node_time, latest_node_time]
    t_min: f64,
    t_max: f64,
}
```

Convert the TBP breakpoints to calendar time once at construction:
`t_calendar = present_time - t_tbp`. Because TBP increases going into the
past while calendar time decreases, the array order reverses; store the
calendar-time arrays in ascending order.

The evaluation method computes the plain-probability coalescent weight at a
single point for a node with a given multiplicity:

```
fn eval(t: f64, multiplicity: f64) -> f64:
    i_t      = integral_merger_rate.eval(t)
    k_t      = lineage_counts.eval(t)
    tc_t     = tc.eval(t)
    lambda_t = k_t * max(k_t - 1, 0.5) / (2 * tc_t)
    neg_log  = multiplicity * (i_t - ln(lambda_t))
    return exp(-neg_log)
```

All arithmetic is scalar; no allocations needed.

### Apply the contribution after combining child messages

Remove the coalescent contribution from the initial seed of `result` in the
backward pass. Instead, apply it as a final pointwise step after all child
messages have been accumulated:

```
// Accumulate child messages as today
result = None
for child in children:
    parent_message = child_time_dist ⊛ (−branch_dist)
    result = result * parent_message   // all Function × Function

// Apply coalescent contribution in-place on the result's grid
if coalescent is Some and result is Some(Function):
    multiplicity = n_children - 1
    result.y *= coalescent.eval_on_grid(result.grid, multiplicity)
    result = result.normalize()
```

`eval_on_grid` evaluates `eval(t, multiplicity)` at each grid point, returning
an array. The multiplication is pointwise on the existing array — no new grid
construction, no overlap computation, no resampling, no `Formula` type
involved at all.

This also removes the dependency on `multiply_formula_function` from the
coalescent path entirely.

### Pass `CoalescentModel` directly into the backward pass

The backward pass signature changes from:

```
propagate_distributions_backward(
    graph,
    coalescent_contributions: Option<&IndexMap<NodeKey, Arc<DistributionNegLog>>>
)
```

to:

```
propagate_distributions_backward(
    graph,
    coalescent: Option<&CoalescentModel>
)
```

The `IndexMap` and all the `Arc<DistributionNegLog>` allocations are
eliminated. The `compute_node_contributions` function and the
`DistributionFormula`-based closure infrastructure can be removed entirely.

### Leaf contributions in the proposed architecture

Leaves have `n_children = 0`, so `multiplicity = max(n_children - 1, 0) = 0`.

The `eval` formula at multiplicity 0 reduces to:

```
neg_log = 0 * (i_t - ln(lambda_t)) = 0
eval(t, 0) = exp(0) = 1
```

This means the `λ(t)` term vanishes at leaves (as expected — leaves are not
coalescent events), but the `exp(I(t))` survival term remains and is encoded
in the `i_t` component. In the refactored architecture the leaf contribution
is handled as follows:

The backward pass for a leaf node produces `result = None` (leaf time
distribution comes from the date constraint, not backward accumulation).
However, the leaf's `exp(-I(t_leaf))` term still needs to flow up to the
parent as a scalar weight. The cleanest approach is:

1. For each leaf, evaluate `coalescent.eval_leaf_scalar(t_leaf)` which returns
   `exp(-I(t_leaf))` (using `i_t` only, no lambda term).
2. Accumulate these scalars into a running log-likelihood total rather than
   incorporating them into the distribution messages.
3. The parent's `exp(I(t_parent))` term is already handled by the parent's
   own `multiplicity = n_children - 1` evaluation in the backward pass.

Alternatively, if tracking a scalar log-likelihood accumulator is inconvenient,
the leaf contribution can be applied to the parent's combined message after
all children are merged, as an additional pointwise factor `exp(-I(t))` with
weight 1 per leaf child. This is equivalent because the leaf's date is fixed
and `exp(-I(t_leaf))` is a constant with respect to the parent's time.

Either way: **do not compute a `DistributionFormula` for leaf nodes**. The
contribution is a scalar evaluated at the leaf's known date, not a function
of the parent's unknown time.

---

## Summary of changes

| Current | Proposed |
|---|---|
| N `DistributionFormula` closures, one per node | One `CoalescentModel` struct |
| TBP→calendar conversion on every eval call | Calendar-time arrays built once at construction |
| `Array1::from_vec` allocations per scalar eval | Scalar arithmetic only |
| Formula applied before child messages, grid guessed | Applied after children, grid is known |
| `multiply_formula_function` path (overlap + resample) | Pointwise multiply on existing array |
| `IndexMap<NodeKey, Arc<DistributionNegLog>>` | Removed entirely |
| Leaf contributions computed but silently discarded (bug) | `eval_leaf_scalar(t_leaf)` accumulated into log-likelihood |
