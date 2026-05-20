# Coalescent backward pass missing leaf and root contributions

v1's backward pass applies only internal node coalescent contributions. Leaf survival contributions and root correction are computed but not consumed (leaves) or not computed at all (root correction). v0 applies all three pieces.

## Problem

The Kingman coalescent neg-log-likelihood decomposes into three per-node pieces via algebraic telescoping of branch survival integrals (see [algo inventory](../algo/timetree.md#kingman-coalescent)):

| Piece           | Neg-log value           | v0                                           | v1                                                                      |
| --------------- | ----------------------- | -------------------------------------------- | ----------------------------------------------------------------------- |
| Internal node   | `m * (I(t) - ln(λ(t)))` | Applied (`clock_tree.py:505-506`)            | Applied (`backward_pass.rs:48-51`)                                      |
| Leaf            | `-I(t_leaf)`            | Applied (`clock_tree.py:474-480`, `499-503`) | Computed but **not applied** (`backward_pass.rs:45`: `None` for leaves) |
| Root correction | `+I(t_root)`            | Applied (`clock_tree.py:518-530`)            | **Not computed**                                                        |

### v0 reference

v0 applies leaf contributions in two code paths during the backward pass:

1. Leaves with precise dates (`clock_tree.py:474-480`): adds `-I(date_peak)` as a constant to the branch length distribution. The constant doesn't change the distribution shape but affects relative weighting between siblings at the parent.

2. Leaves with uncertain dates (`clock_tree.py:499-503`): adds `-I(time_points)` as a distribution to the messages being multiplied. This changes the shape within the date range, biasing the effective date.

3. Root correction (`clock_tree.py:518-530`): after combining all child messages at the root, multiplies by `Distribution(x, +I(x), is_log=True)`. The comment reads: "Removed merger rate must be added back at the root as no longer an internal node."

### v1 current state

`compute_node_contributions()` in [contributions.rs](../../packages/treetime/src/coalescent/contributions.rs) computes leaf and internal contributions correctly (validated by golden master tests). The backward pass at [backward_pass.rs](../../packages/treetime/src/timetree/inference/backward_pass.rs) skips leaf nodes at line 45 (`result = None` for leaves) and does not add root correction.

## Mathematical background

The full coalescent log-likelihood for a tree:

```
ln P(tree | Tc) = Σ_{mergers} m_i * ln(λ(t_i)) - Σ_{branches} (I(t_parent) - I(t_child))
```

The survival sum telescopes: each node j with k_j children contributes `+k_j * I(t_j)` as parent of k_j branches, and `-I(t_j)` as child of its parent's branch (non-root only). Grouping by node:

- Internal (non-root): net coefficient `k_j - 1 = m_j`. Combined with merger: `m*(ln(λ) - I(t))`. Neg-log: `m*(I(t) - ln(λ))`.
- Root: net coefficient `k_root` (no child subtraction). Combined with merger: `m*ln(λ) - k*I(t)`. Equals `m*(ln(λ) - I(t)) - I(t)`. So: internal formula plus root correction `+I(t)`.
- Leaf: net coefficient `-1` (child only, zero children). Neg-log: `-I(t_leaf)`.

The three pieces sum to the exact Kingman neg-log-likelihood.

## Impact

Root correction (most significant): adds `I(t_root)` to root's neg-log, penalizing ancient root times proportionally to the cumulative merger rate. Without it, the coalescent prior does not constrain the root age as much as v0 does. The root time bias propagates to all internal nodes through the forward pass and affects Tc optimization (which depends on node time estimates). Effect magnitude is proportional to `1/Tc`: strong for small Tc (strong coalescent), negligible for large Tc (weak coalescent).

Leaf contributions for precise dates (low impact): `-I(t_leaf)` is a constant that doesn't change the distribution shape. It only affects relative weighting between siblings at the parent node. Second-order effect.

Leaf contributions for uncertain dates (moderate impact): `-I(t)` varies across the date range, biasing the effective date toward more recent times (higher neg-log penalty for ancient dates within the range). Relevant for samples with uncertain collection dates (common in outbreak settings).

Isochronous sampling (no impact): when all tips are sampled at the present, I(0) = 0, so all leaf contributions are zero.

## Proposed solutions

### S1: Apply leaf contributions in backward pass

In `propagate_distributions_backward_single_node()`, when the current node is a leaf with a coalescent contribution, multiply the leaf's `time_distribution` with the contribution (converted to plain) before the leaf distribution is propagated to its parent through convolution. This matches v0's approach at `clock_tree.py:499-503`.

For leaves with precise dates (delta distributions), the multiplication has no effect on shape, only on the peak value. The framework should handle this correctly through its peak normalization.

### S2: Add root correction

After `propagate_distributions_backward()` completes, if coalescent contributions are active, multiply the root's `time_distribution` by a correction factor. The correction is `Distribution::Formula` evaluating `+I(t_tbp)` in neg-log space (or equivalently `exp(-I(t_tbp))` in plain space), with domain matching the root's distribution.

This requires access to `integral_merger_rate` at the backward pass call site. Either pass it as a parameter or compute the correction inside `compute_coalescent_contributions()` and return it alongside the per-node map.

### S3: Validate impact empirically

Before implementing S1 and S2, run v0 and v1 with `--coalescent` on the same heterochronous datasets. Compare root times and internal node times. If the differences are within the convergence tolerance, the missing pieces may be deferred. If root times differ by more than 1% of tree depth, prioritize the root correction (S2).

Datasets to test: flu/h3n2/20 (heterochronous, moderate tree), ebola (heterochronous, short time span), dengue/500 (larger tree).

## Affected commands

- `timetree` with `--coalescent` (directly)
- `timetree` with `--coalescent-opt` (indirectly, Tc optimization depends on node times)
- `timetree` with `--skyline` (indirectly, skyline Tc estimation depends on node times)

## Related issues

- Source: [M-timetree-coalescent-missing-leaf-and-root-contributions.md](../issues/M-timetree-coalescent-missing-leaf-and-root-contributions.md) -- delete after full resolution
- Coalescent total likelihood is emitted now, but it is still incomplete when leaf and root contributions are omitted from the backward pass
- [--coalescent-opt alone skips initial Tc pass](../issues/M-timetree-coalescent-opt-skips-initial.md) -- interaction: initial Tc estimate quality depends on completeness of coalescent prior
- [Positional likelihood metric differs from v0](../issues/M-timetree-positional-likelihood-metric.md) -- both issues affect how convergence metrics reflect the time-tree objective
