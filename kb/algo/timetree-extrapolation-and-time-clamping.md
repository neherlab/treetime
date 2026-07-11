# Timetree: Distribution Extrapolation and Node-Time Monotonicity

This document specifies three related requirements for correct timetree inference: how grid functions should behave outside their support, how intersection boundaries should be computed when combining distributions, and how to enforce the topology constraint that child nodes are younger than their parents.

---

## 1. Extrapolation outside grid support

`GridFn<T>` represents a piecewise-linear function on a finite uniform grid `[x_min, x_max]`. Two policies are needed depending on context:

- **`Zero`** — return `0.0` for any query outside `[x_min, x_max]`. This is the correct default for any bounded probability distribution: leaf date constraints, branch-length likelihoods, and their products are zero outside their stated support.

- **`Constant`** — return the nearest boundary value (`y[0]` to the left, `y[n-1]` to the right). Use this when the distribution is genuinely uninformative beyond the grid edge, i.e. the tail should be treated as flat rather than absent.

Add a `BoundaryBehavior` enum with these two variants and store independent `left_extrap` / `right_extrap` fields on `GridFn`. Default both to `Zero`. Expose builder methods `with_left_extrap` / `with_right_extrap` that return `Self`. Propagate these fields in `mapv`, `resample`, and `negate_arg_inplace` (the last must swap left↔right because negating the argument reflects the domain).

Expose the same builder methods on `DistributionFunction` (delegating to the inner `GridFn`) and on `Distribution<Y>` (no-op for non-Function variants).

### Applying the policies in the inference passes

The backward and forward passes require different tail behaviour on each side:

| Pass                         | Left tail (far past) | Right tail (far future) |
| ---------------------------- | -------------------- | ----------------------- |
| **Backward** (leaves → root) | `Constant`           | `Zero`                  |
| **Forward** (root → leaves)  | `Zero`               | `Constant`              |

**Rationale:**

In the backward pass, each `parent_message` is computed as

```
parent_message = child_time_dist ⊛ (−branch_dist)
```

This message represents "when could the parent be, given this child?" The parent could be arbitrarily far in the past — there is no upper bound on ancestral age imposed by the child alone. The left tail must therefore be `Constant`, not `Zero`. The right tail is `Zero` because the child's sampling date provides a hard upper bound: the parent cannot be more recent than the child.

In the forward pass, `dist_from_parent` is computed as

```
dist_from_parent = parent_except_subtree ⊛ branch_dist
```

This message represents "when could this node be, given its parent and branch?" The node must be after the parent (branch lengths are non-negative) but there is no lower bound from the parent side alone on how far in the future the child could be. The right tail must be `Constant`. The left tail is `Zero` because the parent's time provides a hard lower bound.

Apply the policy immediately after computing each message, before storing or using it further.

---

## 2. Exact intersection boundaries in distribution multiplication

When two distributions are multiplied pointwise, the result is only non-zero on the intersection of their supports. The result grid must span exactly `[max(a.x_min, b.x_min), min(a.x_max, b.x_max)]`, with that intersection treated as empty if `overlap_min >= overlap_max`.

**Do not** filter the existing grid points of either input to those falling inside the overlap, then use the closest existing point as the boundary. This snaps the boundary to the nearest grid point and loses the fractional part of the intersection. When a range boundary falls between two grid points of a function, the snapped boundary produces a result that is either too wide (including a region where one factor is zero) or too narrow (excluding a region where both are non-zero).

Instead, compute `overlap_min` and `overlap_max` from the analytical support boundaries, then resample or evaluate both factors at `n_points` uniformly spaced across the exact `[overlap_min, overlap_max]` interval.

For `Function × Function`, choose `n_points` from the intersection width and the finer of the two grid spacings:

```
dx       = min(a.dx, b.dx)
n_points = round((overlap_max - overlap_min) / dx) + 1
n_points = clamp(n_points, 2, MAX_GRID_POINTS)
```

This replaces the naive `max(a.len, b.len)` choice, which uses the larger of the two input lengths regardless of how much smaller the intersection is than either input's domain. The naive approach inflates grid sizes whenever the two distributions are wide but only partially overlapping.

The same principle applies to `Range × Function` and `Formula × Function`: compute the exact `[overlap_min, overlap_max]` and derive `n_points` from `b.dx` and the overlap width rather than inheriting `b.len`.

For `Function × Function` and `Range × Function` this applies to both `multiply` and `divide`.

---

## 3. Child-time monotonicity clamp

### Why message-passing alone is insufficient

After the forward pass refines each internal node's time distribution and extracts the argmax, the result is not guaranteed to satisfy `child_t >= parent_t`. For near-identical sequences the branch-length distribution is broad and carries little temporal information. The backward message from the subtree (driven by dated leaves) can dominate the combined distribution and pull its peak to a time earlier than the parent's inferred time, producing a negative branch length in the output.

The extrapolation fix above reduces how often this happens — parents are assigned more appropriate (older) times when backward messages are not artificially truncated on the left — but does not eliminate it. When branch lengths genuinely carry no temporal signal, the subtree dates dominate regardless of extrapolation policy. The two fixes address different failure modes and both are necessary.

### Implementation

After extracting the argmax from the combined distribution (`set_likely_time`), apply a monotonicity clamp to every non-root, non-leaf internal node:

```
if child_t < parent_t:
    child_t = parent_t + epsilon      # epsilon ~ 1e-10
```

Leaves are excluded because their times come from date constraints. The root is excluded because it has no parent. A small epsilon offset ensures the child is strictly after the parent, avoiding zero-length branches.

This matches Python TreeTime v0 behaviour (`clock_tree.make_time_tree`).
