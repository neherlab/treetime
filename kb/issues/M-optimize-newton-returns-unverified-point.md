# Newton per-edge optimizer returns points it never verified

> [!IMPORTANT]
> **Decision required.** The Newton method (`--opt-method newton`) can return a branch length with a lower log-likelihood than its start, and it does not report non-convergence. The options below differ in how far the method departs from its current outputs. Established tools safeguard Newton with a derivative bracket, and compare likelihoods only in some modes (see "Comparison with other implementations"). v0 has no Newton method, so v0 supplies no oracle for Newton outputs; options 1 and 2 change Newton outputs, so they need approval. Default Brent outputs stay unchanged.

## Problem

`fn newton_sqrt_inner()` [packages/treetime/src/optimize/method_newton.rs#L14-L49](../../packages/treetime/src/optimize/method_newton.rs#L14-L49) runs Newton-Raphson on $\ell(s^2)$ with $s = \sqrt{t}$, where $t$ is the branch length and $\ell$ the per-edge log-likelihood (substitutions plus the Poisson indel term). It has four exits, and none of them guarantees that the returned point has a log-likelihood at least as high as the start:

- **Non-concave start**: when $\partial^2 \ell / \partial s^2 \geq 0$ at the start, it returns the result of the 100-point grid search `fn grid_search_inner()` [packages/treetime/src/optimize/zero_boundary.rs#L67-L95](../../packages/treetime/src/optimize/zero_boundary.rs#L67-L95) without refining it [method_newton.rs#L27-L29](../../packages/treetime/src/optimize/method_newton.rs#L27-L29). The grid evaluates every point and compares its best positive point with zero, but the start is not a grid point, so the result can be worse than the start
- **Converged**: the step is below `newton_tolerance_sqrt(s)` [method_newton.rs#L34-L36](../../packages/treetime/src/optimize/method_newton.rs#L34-L36), and the loop returns the proposed point, which it never evaluated
- **Curvature break**: at a new point with $\partial^2 \ell / \partial s^2 \geq 0$ the loop breaks and returns that point [method_newton.rs#L40-L45](../../packages/treetime/src/optimize/method_newton.rs#L40-L45). The function evaluated $\ell$ at this point but uses only the derivatives, so it never compares the point with the start
- **Iteration limit**: after `NEWTON_MAX_ITER = 10` iterations it returns the last proposed step, which it never evaluated [method_newton.rs#L48](../../packages/treetime/src/optimize/method_newton.rs#L48)

The steps have no line search and no bracket: a full Newton step can overshoot past a maximum, and only the step clamp to $[s - \sqrt{s^2 + 1},\ s]$ limits it.

## Evidence

Measurements were taken with throwaway instrumentation of `fn optimize_edge()` on `rust` at `77dd3a84`. For every edge and optimize iteration, the probe recorded the exit, the start, the solver output, the final output after `fn reconcile_zero_boundary()`, and a reference maximum: $\ell$ on 800 log-spaced points over $[10^{-4} \cdot \text{one\_mutation}, 10]$ plus $t = 0$, with golden-section refinement around each grid local maximum.

### Example datasets

Measured: `treetime optimize` with default settings on `flu/h3n2/{20,200}`, `ebola/{20,100}`, `zika/86`, `rsv/a/100`, for `--dense true|false` and `--model infer|jc69` (24 Newton runs, 9,544 per-edge optimizations).

| Exit              | Count | Final $\ell$ below reference by more than $10^{-6}$ | Largest shortfall |
| ----------------- | ----: | --------------------------------------------------: | ----------------: |
| Converged         | 8,780 |                                                   0 |   $1\cdot10^{-9}$ |
| Zero shortcut     |   745 |                                                   0 |                 0 |
| Non-concave start |    17 |                                                  15 | $1.5\cdot10^{-3}$ |
| Iteration limit   |     2 |                                                   0 |                 0 |
| Curvature break   |     0 |                                                   - |                 - |

- **No exit returned a point worse than the start** on these datasets
- **The non-concave-start exit loses precision.** The grid has about 30 points per decade and is not refined, so its best point can be up to about 5% away from the maximum (`flu/h3n2/200`: grid $7.28\cdot10^{-4}$, reference $7.11\cdot10^{-4}$)
- **The non-concave-start exit can return zero when a positive maximum exists below the grid.** On 4 edges (`ebola/100`, start $t = 0$) the maximum lies at $t \approx 2.4\cdot10^{-6}$, below the grid lower bound $0.1 \cdot \text{one\_mutation} \approx 5.3\cdot10^{-6}$, so the grid search returns zero with a shortfall of $1.5\cdot10^{-3}$
- **The two iteration-limit exits** started at $t = 4.6\cdot10^{-4}$ and ended at zero after `reconcile_zero_boundary` sent them to the grid search; the final point equals the reference

### Non-unimodal test surface

The K80 surface of [Dinh and Matsen 2017](https://doi.org/10.1214/16-AAP1240) used in `test_dispatch_zero_boundary.rs` (`fn make_dinh_matsen_k80_contribution()`) has a local maximum at $t \approx 0.167$ and its global maximum at $t \approx 6.37$ ($\ell = -2.77254$ against $-2.77333$).

| Start $t_0$ | Newton exit       | Newton output | $\ell$ at output | $\ell$ at start | Final after reconcile |
| ----------: | ----------------- | ------------: | ---------------: | --------------: | --------------------: |
|        0.05 | curvature break   |          1.05 |        -2.774139 |       -2.773405 |                 0.162 |
|         0.6 | curvature break   |             0 |        -2.773498 |       -2.773753 |                 0.163 |
|           2 | non-concave start |          3.01 |        -2.773138 |       -2.773837 |                  3.01 |
|           5 | converged         |          6.37 |        -2.772536 |       -2.772579 |                  6.37 |

- From $t_0 = 0.05$ the curvature-break exit returns $t = 1.05$, **worse than the start**. The final result is correct only because `reconcile_zero_boundary` sees that $\ell(0) > \ell(1.05)$ and runs the grid search. A design that compares the solver output only with the endpoint $t = 0$ would return zero here
- From $t_0 = 2$ the grid returns its own upper end $1.5 t_0 + \text{one\_mutation} = 3.01$, so the true maximum is outside the grid (see [N-optimize-reconcile-grid-extent-unverified.md](N-optimize-reconcile-grid-extent-unverified.md))

### Not measured

> [!IMPORTANT]
> **Investigation required.** The measurement stopped before these runs: `mpox/clade-ii/100` and `dengue/100` (both methods), datasets above 200 tips, models `k80`, `hky85`, and `tn93` with fixed parameters, and the `timetree` command, which calls the same per-edge optimizer. Collect the exit counts and shortfalls of the table above for these settings, in particular whether the curvature-break exit occurs on real data.

## Comparison with other implementations

- **v0** has no Newton branch-length optimizer. `GTR.optimal_t_compressed()` [packages/legacy/treetime/treetime/gtr.py#L816](../../packages/legacy/treetime/treetime/gtr.py#L816) uses only scipy Brent (see [M-optimize-brent-setup-differs-from-v0.md](M-optimize-brent-setup-differs-from-v0.md))
- **IQ-TREE** (`Optimization::minimizeNewton()`, `utils/optimization.cpp`, revision `63c330d9`) keeps a bracket $[x_l, x_h]$ from the sign of the derivative, takes a bisection step when the function is not concave or the Newton step leaves the bracket, and stops after `maxNRStep = 100` steps. `PhyloTree::optimizeOneBranch()` (`tree/phylotree.cpp`) compares likelihoods only when the result exceeds $0.95$ of the maximum branch length: it then restores the original length if the original likelihood is higher. Elsewhere it accepts the Newton result unchecked
- **RAxML-NG** through coraxlib (`corax_opt_minimize_newton_multi()`, `src/corax/optimize/newton.c`, revision `5849dd2`) also keeps a derivative-sign bracket and caps each step at $x_\text{max} / \text{max\_iters}$. It reports `CORAX_OPT_ERROR_NEWTON_LIMIT` at the iteration limit. The caller (`src/corax/optimize/opt_branches.c`) computes the edge log-likelihood of the result and restores the original branch length when the result is worse only in mode `CORAX_OPT_BLO_NEWTON_SAFE` (`check_loglh_improvement()`). RAxML-NG defaults to `CORAX_OPT_BLO_NEWTON_FAST` (`src/Options.cpp`), which skips the check, and switches to the safe mode after a Newton failure (`src/TreeInfo.cpp`)

## Options

1. **Safeguarded Newton with a likelihood check**: keep a bracket from the derivative sign, take bisection steps when the curvature has the wrong sign or the step leaves the bracket, evaluate every returned candidate, and keep the best feasible evaluated point including the start. Report non-convergence (iteration limit) to the caller. Combines the bracket safeguards of IQ-TREE and coraxlib with the likelihood comparison that coraxlib applies in its safe mode. Changes Newton outputs wherever the current search is wrong or imprecise, and a changed search can also select another basin of a multimodal surface
2. **Likelihood check only**: keep the current iteration, evaluate the candidate at every exit, and keep the better of candidate and start. Smallest change; prevents regressions, and the non-concave-start exit keeps its grid precision loss
3. **Remove the Newton method**: Brent is the default and the v0-compatible method; Newton has no v0 counterpart. Removes the defect by removing the CLI value `newton` (a supported method) and its tests

Neither a derivative bracket nor an improvement check proves a global maximum on a multimodal surface: on the Dinh and Matsen K80 surface above, a method that stays near its start finds the local maximum $t \approx 0.167$, not $t \approx 6.37$.

## Cross-references

- [N-optimize-zero-boundary-handled-in-several-places.md](N-optimize-zero-boundary-handled-in-several-places.md): the reconcile step currently hides the curvature-break defect
- [N-optimize-per-edge-constants-without-provenance.md](N-optimize-per-edge-constants-without-provenance.md): `NEWTON_REL_TOL`, `NEWTON_ABS_TOL`, `NEWTON_MAX_ITER`
- [N-optimize-grid-search-resolution-unverified.md](N-optimize-grid-search-resolution-unverified.md): grid precision of the non-concave-start exit
