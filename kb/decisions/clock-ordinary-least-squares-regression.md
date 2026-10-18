# Clock regression is ordinary least squares with an exact root search

The root-to-tip regression of `clock`, `timetree` and the min-dev rooting of `optimize` is ordinary least squares
(OLS) with one unit weight per dated tip, its root search is exact on every edge, and the rate uncertainty comes from
a Gaussian approximation of the timetree, not from the regression.

**Type**: Scientific change (v0 divergence) and removal of options (`--covariation`, `--tip-slack`, the variance and
split-optimizer flags).

**Origin**: phases 0a to 0c of
[timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md),
following `rtt_ols.py` and `gaussian_timetree.py` of
[neherlab/gaussian-timetree](https://github.com/neherlab/gaussian-timetree).

**v1 locations**:

- regression sums and their shift along a branch: `fn ClockSet::shift()`
  [packages/treetime/src/clock/clock_set.rs](../../packages/treetime/src/clock/clock_set.rs)
- backward and forward passes: `fn clock_regression()`
  [packages/treetime/src/clock/clock_regression.rs](../../packages/treetime/src/clock/clock_regression.rs)
- root search: `fn find_best_root()` and `fn optimize_exact()`
  [packages/treetime/src/clock/find_best_root/](../../packages/treetime/src/clock/find_best_root/)
- min-dev rooting of `optimize`: `fn reroot_min_dev()`
  [packages/treetime/src/reroot/orchestrate.rs](../../packages/treetime/src/reroot/orchestrate.rs)
- rate uncertainty: `fn gaussian_clock()` [packages/treetime/src/clock/gaussian.rs](../../packages/treetime/src/clock/gaussian.rs)

## Decision

**Regression.** For a root $r$, the clock fits $d_i(r) = \alpha + \mu t_i$ to the dated tips $i$ that are not
outliers, each with weight one, from the sums $(n, \sum t, \sum t^2, \sum d, \sum d^2, \sum t d)$. Moving the root
by $\ell$ along a branch lengthens the divergences of all tips below it by $\ell$:
$\sum d \mathrel{+}= \ell n$, $\sum td \mathrel{+}= \ell \sum t$, $\sum d^2 \mathrel{+}= 2\ell\sum d + \ell^2 n$.
One backward and one forward pass give the sums on both sides of every edge. Outliers and undated tips contribute
no sums.

**Root search.** With the root at fraction $x$ of an edge of length $\ell$, the score
$\chi^2(x) = \tfrac12 \mathrm{RSS}$ (estimated rate) or $\tfrac{1}{2n}\mathrm{RSS}$ at a fixed rate is a quadratic
in $x$: the date sums do not change and $\sum d$, $\sum td$ are linear and $\sum d^2$ quadratic in $x$. The parabola
through $x = 0, \tfrac12, 1$ is exact, and its vertex clamped to $[0, 1]$ is the optimum on that edge. Every edge is
searched and the minimum over the root and all edges wins. A required positive rate is checked on each edge's
optimum rather than imposed as a constraint, so data without temporal signal stop with the error that asks for a
rate instead of fitting a rate near zero. Clock `min-dev` and the min-dev rooting of `optimize` are the same search
at a fixed rate of zero ([clock-min-dev-fixed-zero-rate.md](clock-min-dev-fixed-zero-rate.md)); `optimize` gives
every tip the same date for it.

**Rate uncertainty.** OLS treats the tips as independent, while they share ancestry, so its covariance is no rate
uncertainty; the unscaled inverse normal matrix that v1 used to report is not one either. The standard deviation of
the rate comes from the Gaussian timetree at the final rate: elapsed time on a branch with $d$ substitutions is
$\mathcal N(\nu/\eta, \nu/\eta^2)$ with $\nu = d + 1$ and $\eta = \mu (D + M + N)/D$, and the observed information
of the likelihood with the root time free gives the rate and root-time variances. It needs the sequence length
(alignment or `--sequence-length`), is reported as `rate_std` of an estimated rate in the clock-model JSON and the
augur node data, and enters `--confidence` unless `--clock-std-dev` is given. A rate the user fixed has no estimated
uncertainty. A timetree rate set by the marginal-likelihood step counts as estimated, with $\chi^2$ and $r$ of the
final clock line.

## v0 divergence

- **No covariation.** v0 `--covariation` weights the regression by branch variances
  $(\ell_e + [e\ \text{terminal}]\, s^2/L)/L$ (tip slack $s$, sequence length $L$). v1 removes the option, and
  `--tip-slack` with it, without a remapping: the shared ancestry that the weights model enters through the Gaussian
  timetree instead, and in the maintainers' experience OLS is the more robust starting point on real data. No
  comparison of the two fits on the validation data was made.
- **Terminal branches.** Without covariation v0 gives each tip variance 1 and splits it along a terminal branch, $x$
  on the tip's side and $1 - x$ on the rest, so a tip's weight grows without bound as the root approaches it. v1
  keeps unit tip weights for every root position, the OLS of `rtt_ols.py`.
- **Root search.** v0 scores a 6-point grid per edge and refines with bounded Brent (`xatol` $10^{-6}$); v1 places
  the root at the exact optimum of every edge.
- **Rate uncertainty.** v0 reports the covariance of the covariation-weighted fit and no rate uncertainty without
  it; v1 reports the Gaussian rate standard deviation and no regression covariance.

## Evidence

- On the 240 validation sets (`treetime clock` with rerooting), the exact search lowers $\chi^2$ in every run, the
  mean $\lvert\log\rvert$ rate error at $L = 10^4$, $\mu = 10^{-2}$ drops from 0.039 to 0.022, and the runs that
  stopped on a negative rate (one without the clock filter, two with it) complete.
- On `data/dengue/100` the fit equals an independent OLS root search (`rtt_ols.py`, rerooting on every node and
  branch) on the tree with the outliers pruned: rate $6.460\cdot10^{-4}$, root 1623.80, $r = 0.9667$. Without the
  filter, v1 before the exact search rejected every root for a negative rate; the exact search finds a positive one
  ($6.88\cdot10^{-4}$).
- The Gaussian rate standard deviation matches the inverse of a finite-difference Hessian of the likelihood. On the
  validation data, $\lvert\hat\mu - \mu\rvert < 1.645\,\sigma_\mu$ in 77-83% of the runs at $L = 10^4$,
  $\mu \ge 10^{-3}$ and in 47% at $L = 10^3$, $\mu = 10^{-2}$ (nominal 90%); at low signal the rate is biased far
  beyond its width (proposal, phase 0b result).

## Validation contract

- On every edge the search returns the analytic optimum of the quadratic score (test: the root-CD edge of a 4-tip
  tree at split $6.18/45$).
- Min-dev rooting of a star $(A{:}1, B{:}2, C{:}6)$ lands $9/4$ up the branch to $C$.
- `ClockSet::shift` equals summing the lengthened divergences directly.
- The dengue/100 golden master matches the independent OLS search.
- `rate_std` matches the inverse finite-difference Hessian of the explicit Gaussian likelihood.
