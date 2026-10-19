# Clock regression is ordinary least squares with an exact root search

The root-to-tip regression of `clock`, `timetree` and the min-dev rooting of `optimize` is ordinary least squares
(OLS) with one unit weight per dated tip, its root search is exact on every edge, and the rate uncertainty comes from
a Gaussian approximation of the timetree, not from the regression.

**Type**: Scientific change (v0 divergence) and removal of options (`--covariation`, `--tip-slack`, the variance and
split-optimizer flags).

**Origin**: phases 0a to 0c of
[timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md),
following an unpublished derivation and reference scripts of the maintainers (`rtt_ols.py` for the regression,
`gaussian_timetree.py` for the Gaussian timetree). Section "Gaussian timetree" below states the formulas that v1 uses.

**v1 locations**:

- regression sums and their shift along a branch: `fn ClockSet::shift()`
  [packages/treetime/src/clock/clock_set.rs](../../packages/treetime/src/clock/clock_set.rs)
- backward and forward passes: `fn clock_regression()`
  [packages/treetime/src/clock/clock_regression.rs](../../packages/treetime/src/clock/clock_regression.rs)
- root search: `fn find_best_root()` and `fn optimize_exact()`
  [packages/treetime/src/clock/find_best_root/](../../packages/treetime/src/clock/find_best_root/)
- reroot driver of every command, and its date-free entry for `optimize`: `fn reroot_at_best_root()`,
  `fn reroot_undated()` [packages/treetime/src/clock/reroot.rs](../../packages/treetime/src/clock/reroot.rs)
- Gaussian timetree: `struct GaussianTree` [packages/treetime/src/clock/gaussian.rs#L107-L121](../../packages/treetime/src/clock/gaussian.rs#L107-L121), built by `fn GaussianTree::new()`
  [packages/treetime/src/clock/gaussian.rs#L124-L196](../../packages/treetime/src/clock/gaussian.rs#L124-L196); node marginals and rate uncertainty `fn GaussianTree::clock()` [packages/treetime/src/clock/gaussian.rs#L198-L226](../../packages/treetime/src/clock/gaussian.rs#L198-L226);
  closed-form rate `fn GaussianTree::rate()` [packages/treetime/src/clock/gaussian.rs#L228-L235](../../packages/treetime/src/clock/gaussian.rs#L228-L235); callers `fn timetree_gaussian_clock()`
  [packages/treetime/src/timetree/init.rs#L152-L180](../../packages/treetime/src/timetree/init.rs#L152-L180) and `fn gaussian_start_rate()` [packages/treetime/src/timetree/init.rs#L184-L200](../../packages/treetime/src/timetree/init.rs#L184-L200)

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
$\mathcal N(\nu/\eta, \nu/\eta^2)$ with $\nu = d + 1$ and $\eta = \mu L (D + M + N)/D$, and the observed information
of the likelihood with the root time free gives the rate and root-time variances (section "Gaussian timetree"). It needs the sequence length
(alignment or `--sequence-length`), is reported as `rate_std` of an estimated rate in the clock-model JSON and the
augur node data, and enters `--confidence` unless `--clock-std-dev` is given. A rate the user fixed has no estimated
uncertainty. A timetree rate set by the marginal-likelihood step counts as estimated, with $\chi^2$ and $r$ of the
final clock line.

## Gaussian timetree

The Gaussian timetree approximates the time inference of TreeTime <a id="cite-1"></a>[Sagulenko et al. 2018](https://doi.org/10.1093/ve/vex042) [[1](#ref-1)]
at a given clock rate by Gaussian <a id="gloss-use-1"></a>messages <sup>[1](#gloss-1)</sup>. v1 uses it for the rate and
root-time uncertainty, for the internal-node marginals that start the timetree loop, and for a start rate when the
root-to-tip regression gives no positive rate. The formulas follow the unpublished derivation of the maintainers.
Like the least-squares dating of LSD <a id="cite-2"></a>[To et al. 2016](https://doi.org/10.1093/sysbio/syv068) [[2](#ref-2)], it is a normal
approximation of the Poisson molecular clock of <a id="cite-3"></a>[Langley and Fitch 1974](https://doi.org/10.1007/BF01797451) [[3](#ref-3)]. LSD
puts the normal error on the branch length, with variance $(b + c/L)/L$ for a branch length $b$ and a smoothing
constant $c$; the Gaussian timetree puts it on the elapsed time of the branch, and its messages give node marginals
and a likelihood of the rate besides the point estimate.

**Branches.** A branch $e$ of length $\ell_e$ carries

$$\nu_e = \ell_e L + 1$$

where:

- $L$ -- sequence length (alignment length or `--sequence-length`)
- $\ell_e L$ -- expected number of substitutions on the branch

The elapsed time on the branch is $\mathcal N(\nu_e/\eta,\ \nu_e/\eta^2)$, the
<a id="gloss-use-2"></a>moment match <sup>[2](#gloss-2)</sup> of its $\mathrm{Gamma}(\nu_e, \eta)$ likelihood under a flat prior on
time, where $\eta$ is the rate of the event process. The $+1$ in $\nu_e$ is the event that ends a branch: a split,
or the sampling of a tip ([packages/treetime/src/clock/gaussian.rs#L144-L152](../../packages/treetime/src/clock/gaussian.rs#L144-L152)).

**Event rate and substitution rate.** The event process counts substitutions and branch ends, so

$$\eta = \mu L \frac{D + M + N}{D}, \qquad \mu = \eta \frac{D}{L (D + M + N)}$$

where:

- $\mu$ -- substitution rate per site and unit time
- $D = \sum_e (\nu_e - 1)$ -- substitutions on the branches that carry a message
- $M$, $N$ -- number of internal and terminal branches that carry a message

A branch carries a message when a dated tip that is not a clock-filter outlier lies below it
(`fn GaussianTree::eta()` and `fn GaussianTree::mu()` [packages/treetime/src/clock/gaussian.rs#L237-L243](../../packages/treetime/src/clock/gaussian.rs#L237-L243)).

**Messages.** A message about the time $t$ of a node is a triple $(S, R, W)$ that does not depend on $\eta$, with

$$\mathbb E[t] = S - \frac{R}{\eta}, \qquad \operatorname{var}(t) = \frac{W}{\eta^2}$$

A dated tip sends $(t_i, 0, 0)$, its date $t_i$ known exactly. A message crosses branch $e$ toward the parent as
$(S, R + \nu_e, W + \nu_e)$ and toward the child as $(S, R - \nu_e, W + \nu_e)$ (`fn Message::toward_parent()`,
`fn Message::toward_child()` [packages/treetime/src/clock/gaussian.rs#L34-L58](../../packages/treetime/src/clock/gaussian.rs#L34-L58)).

**Precision-weighted combination.** Messages $k = 1, \dots, K$ about the same time combine by summing their
<a id="gloss-use-3"></a>natural parameters <sup>[3](#gloss-3)</sup> $(1/W_k, S_k/W_k, R_k/W_k)$:

$$\frac{1}{W} = \sum_k \frac{1}{W_k}, \qquad S = W \sum_k \frac{S_k}{W_k}, \qquad R = W \sum_k \frac{R_k}{W_k}$$

One pass toward the root and one toward the tips give every internal node $v$ the combination $(S_v, R_v, W_v)$ of
all its messages, and its marginal time $\mathcal N(S_v - R_v/\eta,\ W_v/\eta^2)$ (`struct Natural`,
`fn combine()` [packages/treetime/src/clock/gaussian.rs#L61-L99](../../packages/treetime/src/clock/gaussian.rs#L61-L99); `fn backward()`, `fn forward()` [packages/treetime/src/clock/gaussian.rs#L246-L312](../../packages/treetime/src/clock/gaussian.rs#L246-L312)).

**Likelihood of the rate.** Each combination at a node leaves $-\tfrac12 (a\eta^2 - 2b\eta + c)$ in the
log-likelihood, with

$$a = \sum_k \frac{S_k^2}{W_k} - \frac{S^2}{W}, \qquad b = \sum_k \frac{S_k R_k}{W_k} - \frac{S R}{W}, \qquad c = \sum_k \frac{R_k^2}{W_k} - \frac{R^2}{W}$$

With the root time at its optimum, the log-likelihood of the event rate is

$$\log \mathcal L(\eta) = \kappa \ln \eta - \tfrac12 A \eta^2 + B \eta - \tfrac12 C$$

where:

- $\kappa$ -- number of dated tips that carry a message
- $A$, $B$, $C$ -- sums of $a$, $b$, $c$ over the combinations at all internal nodes

$C$ does not depend on $\eta$, so v1 does not compute it.

**Closed-form rate.** $\partial_\eta \log \mathcal L = 0$ is the quadratic $A\eta^2 - B\eta - \kappa = 0$, whose
positive root is

$$\eta^* = \frac{B + \sqrt{B^2 + 4 A \kappa}}{2A}$$

It exists whenever $A > 0$, that is, whenever the tip dates vary along the tree, so it is positive where the
root-to-tip regression may give a negative slope. $\mu^*$ follows from $\eta^*$ by the conversion above
(`fn GaussianTree::rate()` [packages/treetime/src/clock/gaussian.rs#L228-L235](../../packages/treetime/src/clock/gaussian.rs#L228-L235)).

**Rate and root uncertainty.** With the root time $t_r$ fixed, the log-likelihood is
$\kappa \ln \eta - \tfrac12 A' \eta^2 + B' \eta - \tfrac12 C'$ with $A' = A + (t_r - S_r)^2/W_r$, where
$(S_r, R_r, W_r)$ is the combination at the root. Its <a id="gloss-use-4"></a>observed information <sup>[4](#gloss-4)</sup> at
$t_r = S_r - R_r/\eta$ is

$$I = \begin{pmatrix} \kappa/\eta^2 + A + R_r^2/(\eta^2 W_r) & -R_r/W_r \\ -R_r/W_r & \eta^2/W_r \end{pmatrix}$$

with rows and columns $(\eta, t_r)$, and

$$\sigma_\eta = \sqrt{\frac{I_{22}}{\det I}}, \qquad \sigma_{t_r} = \sqrt{\frac{I_{11}}{\det I}}, \qquad \sigma_\mu = \sigma_\eta \frac{D}{L (D + M + N)}$$

Both are undefined unless $\det I > 0$ and $I_{11} > 0$. `rate_std` is $\sigma_\mu$ at the final clock rate
(`fn rate_and_root_std()` [packages/treetime/src/clock/gaussian.rs#L314-L320](../../packages/treetime/src/clock/gaussian.rs#L314-L320)).

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

## Glossary

1. <a id="gloss-1"></a> **Message.** A summary of the evidence about the time of a node that one part of the tree
   sends to another, here a Gaussian given by the triple $(S, R, W)$. [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Moment match.** The Gaussian with the same mean and variance as a given distribution.
   [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **Natural parameters.** The parameters in which the logarithm of a Gaussian density is linear,
   precision and precision times mean; the product of Gaussian densities about the same variable sums them.
   [↩](#gloss-use-3)
4. <a id="gloss-4"></a> **Observed information.** The negative Hessian of the log-likelihood at the estimate; its
   inverse approximates the covariance of the estimated parameters. [↩](#gloss-use-4)

## References

1. <a id="ref-1"></a> Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-likelihood phylodynamic analysis." _Virus Evolution_ 4 (1): vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-1)
2. <a id="ref-2"></a> To, Thu-Hien, Matthieu Jung, Samantha Lycett, and Olivier Gascuel. 2016. "Fast dating using least-squares criteria and algorithms." _Systematic Biology_ 65:82-97. https://doi.org/10.1093/sysbio/syv068 [↩](#cite-2)
3. <a id="ref-3"></a> Langley, Charles H., and Walter M. Fitch. 1974. "An examination of the constancy of the rate of molecular evolution." _Journal of Molecular Evolution_ 3:161-177. https://doi.org/10.1007/BF01797451 [↩](#cite-3)
