# Coalescent backward pass missing leaf and root contributions

v1's backward pass applies only internal node coalescent contributions. Leaf survival contributions and root correction are computed but not consumed (leaves) or not computed at all (root correction). v0 applies all three pieces.

## Problem

On a bifurcating tree, the ordinary Kingman coalescent neg-log-likelihood decomposes into three per-node pieces via algebraic telescoping of branch survival integrals (see [algo inventory](../algo/timetree.md#kingman-coalescent)):

| Piece           | Neg-log value           | v0                                           | v1                                                                      |
| --------------- | ----------------------- | -------------------------------------------- | ----------------------------------------------------------------------- |
| Internal node   | `m * (I(t) - ln(λ(t)))` | Applied (`clock_tree.py:505-506`)            | Applied (`backward_pass.rs:48-51`)                                      |
| Leaf            | `-I(t_leaf)`            | Applied (`clock_tree.py:474-480`, `499-503`) | Computed but **not applied** (`backward_pass.rs:45`: `None` for leaves) |
| Root correction | `+I(t_root)`            | Applied (`clock_tree.py:518-530`)            | **Not computed**                                                        |

### v0 reference

v0 applies leaf contributions in two code paths during the backward pass:

1. Leaves with precise dates (`clock_tree.py:474-480`): adds $-I(t_\mathrm{leaf})$ as a constant to the branch length distribution. It cannot change that leaf message's shape or the relative weighting of sibling times at a fixed parameter value. Retaining the constant affects the accumulated objective and comparisons across $T_c$ values.

2. Leaves with uncertain dates (`clock_tree.py:499-503`): adds $-I(t)$ as a distribution to the messages being multiplied. Because $I(t)$ is nondecreasing in time before present, $-I(t)$ is smaller for older dates. Equivalently, the corresponding likelihood factor is $\exp[I(t)]$, which favors an older child time closer to its parent within the allowed date range.

3. Root correction (`clock_tree.py:518-530`): after combining all child messages at the root, multiplies by `Distribution(x, +I(x), is_log=True)`. The comment reads: "Removed merger rate must be added back at the root as no longer an internal node."

### v1 current state

`compute_node_contributions()` in [contributions.rs](../../packages/treetime/src/coalescent/contributions.rs) computes the leaf and internal terms of the telescoped objective. The backward pass at [backward_pass.rs](../../packages/treetime/src/timetree/inference/backward_pass.rs) skips leaf nodes at line 45 (`result = None` for leaves) and does not add the root correction. `DistributionNegLog` permits finite negative ordinates, but the API does not distinguish a normalized distribution from a coupled additive objective term whose meaning depends on leaf, internal, and root terms being consumed together.

## Mathematical background

For a bifurcating tree, let $t_i$ be merger times, $t_p$ and $t_c$ be parent and child times, $\lambda(t)$ be the total merger rate, $\kappa(t)$ be the per-lineage merger rate, and $I(t)=\int_0^t\kappa(u)\,du$ be its cumulative integral. The full coalescent log-likelihood is

$$
\log P(\mathcal T\mid T_c)
= \sum_{i\in\mathrm{mergers}} \log \lambda(t_i)
- \sum_{(p,c)\in\mathrm{branches}} \left[I(t_p)-I(t_c)\right],
$$

where $\mathcal T$ is the time tree and $T_c$ parameterizes the effective population-size scale.

The survival sum telescopes: each node j with k_j children contributes `+k_j * I(t_j)` as parent of k_j branches, and `-I(t_j)` as child of its parent's branch (non-root only). Grouping by node:

- Internal non-root node: two parent occurrences and one child occurrence give net coefficient one. Combined with its merger event, the neg-log term is $I(t)-\log\lambda(t)$.
- Root: two parent occurrences and no child occurrence give the internal term plus the correction $+I(t_\mathrm{root})$.
- Leaf: one child occurrence and no parent occurrence give the neg-log correction $-I(t_\mathrm{leaf})$.

These three pieces sum to the ordinary Kingman neg-log-likelihood for a bifurcating tree. v0 generalizes an internal node with $k$ children to multiplicity $m=k-1$. Because ordinary Kingman merger events are binary, that multifurcation rule is v0 parity behavior until a binary-resolution or multiple-merger model is approved.

## Impact

**Root correction** (most significant): adds `I(t_root)` to root's neg-log, penalizing ancient root times proportionally to the cumulative merger rate. Without it, the coalescent prior does not constrain the root age as much as v0 does. The root time bias propagates to all internal nodes through the forward pass and affects Tc optimization (which depends on node time estimates). Effect magnitude is proportional to `1/Tc`: strong for small Tc (strong coalescent), negligible for large Tc (weak coalescent).

**Leaf contributions for precise dates**: $-I(t_\mathrm{leaf})$ is a constant that does not change the distribution shape or sibling-time weighting. Omitting it changes the accumulated objective and can affect comparisons across $T_c$ values.

**Leaf contributions for uncertain dates**: $-I(t)$ varies across the date range and favors older child times closer to the parent. Omitting it changes the inferred effective date within the permitted range.

**Isochronous sampling** (no impact): when all tips are sampled at the present, I(0) = 0, so all leaf contributions are zero.

## Decisions required

- Approve a representation for the complete telescoped objective, including signed leaf corrections and the coupled root correction.
- Decide whether contributions are applied atomically during message passing or evaluated as a separate complete coalescent objective.
- Establish analytical tests on bifurcating trees and v0 golden masters for precise-date leaves, uncertain-date leaves, and the root.
- Decide whether multifurcations preserve v0's multiplicity rule, are resolved to binary events, or use an approved multiple-merger model.

No implementation ticket is ready until the contribution API contract is approved.

## Affected commands

- `timetree` with `--coalescent` (directly)
- `timetree` with `--coalescent-opt` (indirectly, Tc optimization depends on node times)
- `timetree` with `--skyline` (indirectly, skyline Tc estimation depends on node times)

## Related issues

- Coalescent total likelihood is emitted now, but it is still incomplete when leaf and root contributions are omitted from the backward pass
- [Positional likelihood metric differs from v0](M-timetree-positional-likelihood-metric.md) - both issues affect how convergence metrics reflect the time-tree objective
