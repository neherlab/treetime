# Input-mode Poisson branch likelihoods diverge from v0 after the first cell and beyond the grid

## Problem

With `--branch-length-mode=input`, `fn InputBranchLikelihood::distribution()` ([packages/treetime/src/timetree/inference/input_branch_length.rs](../../packages/treetime/src/timetree/inference/input_branch_length.rs)) samples v0's branch-length likelihood (`branch_len_interpolator.py`, the `branch_length_mode == 'input'` branch) on a uniform grid with a step of a tenth of the peak width. The first cell of a Poisson branch follows v0's cost exactly, and the Gaussian cost stays within the bound of linear interpolation, $h^2/(8\sigma^2) = 1.25 \cdot 10^{-3}$ nats ([kb/decisions/timetree-input-branch-likelihood-first-cell-and-bounds.md](../decisions/timetree-input-branch-likelihood-first-cell-and-bounds.md)). Two differences of the Poisson cost remain.

$k$ is the number of mutations on the branch, $L$ the sequence length, $c = 10^{-3}/L$ the v0 offset `MIN_BRANCH_LENGTH * one_mutation`, and $h$ the grid step.

### Cells after the first

Between grid points, the negative log likelihood is linear. The term $k \ln(x + c)$ has the curvature $k/(x + c)^2$, at most $k/(ih)^2$ in cell $i \ge 1$, so the error is at most $k/(8i^2)$ nats in cell $i$, independent of $h$. Largest differences at the v0 grid points within 100 nats of the peak, all in cell 1:

| Case | Difference | Bound $k/8$ |
| --- | --- | --- |
| $k = 0.05$, $L = 1000$ | $2.9 \cdot 10^{-3}$ | $6.3 \cdot 10^{-3}$ |
| $k = 1$, $L = 1000$ | 0.059 | 0.125 |
| $k = 3$, $L = 1000$ | 0.138 | 0.375 |
| $k = 10$, $L = 1000$ | 0.587 | 1.25 |
| $k = 30$, $L = 30000$ | 1.76 | 3.75 |

The golden-master test `test_gm_input_branch_length_poisson_differs_from_v0_within_the_bound_of_each_cell` in `packages/treetime/src/timetree/inference/__tests__/test_gm_input_branch_length.rs` asserts the bound at every v0 grid point inside the grid.

Fix direction: a grid that is non-uniform near zero branch length, as v0's, removes the error, because the curvature falls with the distance from 0. A finer uniform grid does not: its cell 1 lies closer to 0, and the bound $k/8$ stays. A non-uniform grid changes the uniform `GridFn` that convolutions, products and the mass profile use ([kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md)).

### Beyond the grid

The Poisson grid ends where the cost is 100 nats above its minimum (`POISSON_GRID_NEG_LOG_SPAN`) and continues with a linear tail of the slope there. The cost is convex, so the tail underestimates v0 by $k(r - 1 - \ln r)$, with $r$ the branch length divided by the grid end. At v0's last grid point ($x \approx 4$) the relative difference is 0.04% ($k = 0.05$), 0.8% ($k = 1$), 6.4% ($k = 10$), 16% ($k = 30$, $L = 30000$) and 54% ($k = 499$). All these points lie more than 3900 nats above the peak. v0's own Gaussian cost turns linear at 100 nats, and an exact tail to v0's `MAX_BRANCH_LENGTH = 4` needs 40 000 to 220 000 grid points per branch.

## Validation

- Cells after the first: until the grid changes, `test_gm_input_branch_length_poisson_differs_from_v0_within_the_bound_of_each_cell` asserts the accepted bound $k/(8i^2)$ nats in cell $i \ge 1$ at every v0 grid point inside the grid ([kb/decisions/timetree-input-branch-likelihood-first-cell-and-bounds.md](../decisions/timetree-input-branch-likelihood-first-cell-and-bounds.md)). A non-uniform grid resolves this part when every v0 grid point within 100 nats of the peak lies within the interpolation bound $h_i^2 \max |y''| / 8$ of the cell $i$ that holds it, and the test asserts that bound in place of $k/(8i^2)$
- Beyond the grid: the linear tail at 100 nats above the peak is the decided representation (section "Grid reach of the Poisson cost" of the same decision). A tail that follows the convex cost needs a new decision; it resolves this part when the difference at v0's grid points beyond the grid end falls below the bound of the cells inside it
