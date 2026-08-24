# Timetree: Distribution Extrapolation and Committed Node Times

This document specifies how grid functions behave outside their support and how intersection boundaries are computed when combining distributions. It also records the current projection of committed internal-node times; the statistical contract for that projection is unresolved.

---

## 1. Extrapolation outside grid support

`GridFn<T>` represents a piecewise-linear function on a finite uniform grid `[x_min, x_max]`. It stores only the grid and ordinates. Generic interpolation outside the grid returns an error unless the caller supplies explicit left and right `BoundaryBehavior` values.

`DistributionFunction` stores independent `left_extrap` and `right_extrap` fields. Fresh function distributions use `Error` on both sides. `Distribution` exposes both policies through variant-independent optional getters and fits a `Linear` tail from function values near a selected edge. A `Formula` has no stored boundary policy because its exact closure evaluates on demand. `Empty`, `Point`, and `Range` distributions have exact compact support, so their getters return `Hard`, attempts to replace those boundaries have no effect, and evaluation outside their support returns the axis-specific zero-probability value.

The distribution layer supports these policies depending on context:

- **`Hard`** -- return `0.0` for any query outside `[x_min, x_max]`. This is the correct default for any bounded probability distribution: leaf date constraints, branch-length likelihoods, and their products are zero outside their stated support.

- **`Constant`** -- return the nearest boundary value (`y[0]` to the left, `y[n-1]` to the right). Use this when the distribution is genuinely uninformative beyond the grid edge, i.e. the tail should be treated as flat rather than absent.

`DistributionFunction` preserves its policies during resampling and vertical shifts. Scaling transforms fitted-law coefficients, and argument negation reflects each law and swaps the two sides. Formula arithmetic uses its sampling range to construct finite result grids while direct evaluation uses the exact closure.

### Applying the policies in the inference passes

The backward and forward passes require different tail behaviour on each side:

| Pass                          | Left tail (far past) | Right tail (far future) |
| ----------------------------- | -------------------- | ----------------------- |
| **Backward** (leaves -> root) | `Constant`           | `Hard`                  |
| **Forward** (root -> leaves)  | `Hard`               | `Constant`              |

**Rationale:**

In the backward pass, each `parent_message` is computed as

```
parent_message = child_time_dist ⊛ (-branch_dist)
```

This message represents "when could the parent be, given this child?" The parent could be arbitrarily far in the past -- there is no upper bound on ancestral age imposed by the child alone. The left tail must therefore be `Constant`, not `Hard`. The right tail is `Hard` because the child's sampling date provides a hard upper bound: the parent cannot be more recent than the child.

In the forward pass, `dist_from_parent` is computed as

```
dist_from_parent = parent_except_subtree ⊛ branch_dist
```

This message represents "when could this node be, given its parent and branch?" The node must be after the parent (branch lengths are non-negative) but there is no lower bound from the parent side alone on how far in the future the child could be. The right tail must be `Constant`. The left tail is `Hard` because the parent's time provides a hard lower bound.

Apply the policy immediately after computing each message, before storing or using it further.

---

## 2. Exact intersection boundaries in distribution arithmetic

Range-Function, Function-Function, and Formula-Function multiplication are defined on the exact intersection of their operand supports. Range-Function and Function-Function division use the same intersection while the divisor has the default `Error` boundary behavior. A positive-width result grid spans exactly `[max(a.x_min, b.x_min), min(a.x_max, b.x_max)]`. An intersection is empty when `overlap_min > overlap_max`; exact endpoint contact produces a Point distribution evaluated at the shared endpoint, matching v0.

An explicit divisor tail changes the divisor's evaluable domain on that side. `Hard` and `Constant` tails extend division to the dividend boundary; an `Error` side remains constrained by the nominal divisor grid. Left and right sides are resolved independently before intersecting the resulting divisor domain with the dividend. [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md) defines these boundary behaviors and their use by inference messages.

**Do not** filter the existing grid points of either input to those falling inside the overlap, then use the closest existing point as the boundary. This snaps the boundary to the nearest grid point and loses the fractional part of the intersection. When a range boundary falls between two grid points of a function, the snapped boundary produces a result that is either too wide (including a region where one factor is zero) or too narrow (excluding a region where both are non-zero).

Instead, compute `overlap_min` and `overlap_max` from the analytical or explicitly extended boundaries, then resample or evaluate both factors at `n_points` uniformly spaced across the exact `[overlap_min, overlap_max]` interval. Range-Function multiplication and division and Formula-Function multiplication derive `n_points` from the Function spacing. Function-Function multiplication and division use the finer input spacing. Endpoint-inclusive uniform spacing preserves both exact boundaries. [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md) records this intentional divergence from v0's union of knots.

For spacing-derived grids, choose `n_points` from the intersection width and the relevant input spacing:

```
dx       = min(a.dx, b.dx)  # Function-Function; use function.dx for Range-Function or Formula-Function
n_points = clamp(round((overlap_max - overlap_min) / dx) + 1, 2, 1_000_000)
```

Both multiplication and division honor operand tails when computing the support intersection. A `Constant` tail extends the evaluable domain on that side to the other operand's grid boundary. `Hard` and `Error` tails keep the grid boundary as-is (intersection is correct when the value outside support is zero or undefined). This keeps generic `GridFn` evaluation erroring by default and makes every extrapolated arithmetic path explicit. See [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md) for the per-side extension rules.

---

## 3. Committed node-time projection

### Why message-passing alone is insufficient

After the forward pass refines each internal node's time distribution and extracts the argmax, the result is not guaranteed to satisfy `child_t >= parent_t`. For near-identical sequences the branch-length distribution is broad and carries little temporal information. The backward message from the subtree (driven by dated leaves) can dominate the combined distribution and pull its peak to a time earlier than the parent's inferred time, producing a negative branch length in the output.

The extrapolation policy above can change where posterior modes occur, but it does not jointly constrain independently selected modes. Current v1 separately projects committed internal-node point estimates after inference; that projection has a distinct, unresolved statistical contract.

### Current implementation

After extracting the mode from the combined distribution, `set_likely_time()` commits the following value for every non-root, non-leaf internal node:

```
child_t = max(marginal_mode, parent_t)
```

Leaves retain their observed date constraints, and the root has no parent. An inverted internal-node mode is projected to equality with its parent, so the implementation permits zero-duration edges.

This differs from v0 marginal inference, which commits independent posterior modes and can produce negative branch lengths. The projection changes the committed point estimate without recomputing the posterior or its summaries. Its desired contract remains undecided in [kb/issues/M-timetree-marginal-node-times-can-violate-topology.md](../issues/M-timetree-marginal-node-times-can-violate-topology.md); this algorithm inventory records current behavior without approving the divergence.
