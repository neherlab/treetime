# min_dev ignores the fixed slope when scoring candidate roots

## Summary

TreeTime v0 documents `min_dev` as minimizing the variation of <a id="gloss-use-1"></a>root-to-tip distance <sup>[1](#gloss-1)</sup> and passes a fixed zero clock rate into its root search. The returned regression reports that zero rate, but `base_regression()` still scores every candidate with the residual from an estimated-rate temporal regression. The selected root therefore implements the `least-squares` objective rather than the documented `min_dev` objective whenever sampling dates vary.

TreeTime v1 corrects this defect: clock `min-dev` uses `RootObjective::FixedRate(0.0)`, while `least-squares` retains estimated-rate temporal regression.

## Intended objective

For a candidate root under diagonal no-covariation weighting, the <a id="gloss-use-2"></a>weighted least-squares (WLS) <sup>[2](#gloss-2)</sup> objective at a fixed clock rate is

$$
J(\alpha \mid \mu) = \frac{1}{2W}\sum_{i=1}^{n} w_i\left(d_i - \alpha - \mu t_i\right)^2.
$$

Here, $J$ is the normalized residual objective, $\alpha$ is the fitted intercept, $\mu$ is the fixed clock rate, $W = \sum_i w_i$ is total precision, $n$ is the number of contributing tips, $w_i$ is tip $i$'s precision, $d_i$ is its root-to-tip distance, and $t_i$ is its sampling date.

At $\mu = 0$, minimizing over $\alpha$ gives

$$
\widehat{\alpha}_0 = \frac{S_d}{W},
\qquad
J_{\min}(0) = \frac{S_{dd}W - S_d^2}{2W^2}.
$$

Here, $S_d = \sum_i w_i d_i$ is the precision-weighted distance sum and $S_{dd} = \sum_i w_i d_i^2$ is the precision-weighted squared-distance sum. Thus fixed-zero-rate WLS minimizes weighted root-to-tip variance. This is the contract stated by v0 and the minimum-variance rooting criterion described by <a id="cite-1a"></a>[Mai et al. 2017](https://doi.org/10.1371/journal.pone.0182238) [[1](#ref-1)].

## Actual v0 score

`def TreeTime.reroot()` sets `slope = 0.0` for `min_dev` [packages/legacy/treetime/treetime/treetime.py#L575](../../packages/legacy/treetime/treetime/treetime.py#L575), and `def TreeRegression.optimal_reroot()` says that this fixed slope corresponds to minimum root-to-tip variance [packages/legacy/treetime/treetime/treeregression.py#L416-L435](../../packages/legacy/treetime/treetime/treeregression.py#L416-L435).

However, `def base_regression()` uses the supplied slope only when calculating the intercept [packages/legacy/treetime/treetime/treeregression.py#L26-L37](../../packages/legacy/treetime/treetime/treeregression.py#L26-L37). Whenever dates vary, its candidate score remains

$$
J_{\mathrm{v0}} = \frac{1}{2}\left[
S_{dd} - \frac{S_d^2}{W}
- \frac{\left(S_{dt} - S_dS_t/W\right)^2}{S_{tt} - S_t^2/W}
\right].
$$

Here, $J_{\mathrm{v0}}$ is the score returned by `base_regression()`, $S_t = \sum_i w_i t_i$ is the weighted date sum, $S_{tt} = \sum_i w_i t_i^2$ is the weighted squared-date sum, and $S_{dt} = \sum_i w_i d_i t_i$ is the weighted distance-date cross-product. The final subtraction is the reduction in residual obtained by estimating the clock rate. It is independent of the supplied `slope` argument [packages/legacy/treetime/treetime/treeregression.py#L38-L48](../../packages/legacy/treetime/treetime/treeregression.py#L38-L48).

Consequently, `slope=0` changes the reported slope and intercept but does not change candidate ranking from estimated-rate least squares. The command-line help independently states that `min_dev` minimizes root-to-tip variance [packages/legacy/treetime/treetime/argument_parser.py#L78-L87](../../packages/legacy/treetime/treetime/argument_parser.py#L78-L87).

## Scientific distinction

For <a id="gloss-use-3"></a>heterochronous samples <sup>[3](#gloss-3)</sup>, a <a id="gloss-use-4"></a>strict molecular clock <sup>[4](#gloss-4)</sup> predicts

$$
\operatorname{E}[d_i] = \mu\left(t_i - t_r\right),
$$

where $\operatorname{E}[d_i]$ is the expected root-to-tip distance and $t_r$ is the root date. Temporal rooting therefore estimates $\mu$ and minimizes deviations from this line. TempEst describes its best-fit root as the position minimizing squared residuals from the root-to-tip regression <a id="cite-2"></a>[Rambaut et al. 2016](https://doi.org/10.1093/ve/vew007) [[2](#ref-2)], and TreeTime describes the same temporal-rooting criterion <a id="cite-3"></a>[Sagulenko et al. 2018](https://doi.org/10.1093/ve/vex042) [[3](#ref-3)].

Minimum-variance rooting answers a different question: it ignores sampling-date values and seeks the root with the most uniform tip depths. It is justified as a branch-length heuristic when departures from a strict clock are random <a id="cite-1b"></a>[Mai et al. 2017](https://doi.org/10.1371/journal.pone.0182238) [[1](#ref-1)]. For heterochronous data, uniform tip depths are not themselves a strict-clock prediction because later samples have had more time to accumulate substitutions.

## Evidence standard

The erratum is supported independently by:

- the `min_dev` caller supplying `slope=0`;
- two v0 docstrings defining that value as minimum root-to-tip variance;
- `base_regression()` omitting the supplied slope from its score;
- published methods distinguishing temporal-regression rooting from minimum-variance rooting.

## v1 status

Corrected. `fn select_root()` dispatches `RerootMethod::MinDev` to `RootObjective::FixedRate(0.0)` [packages/treetime/src/clock/reroot.rs#L152-L184](../../packages/treetime/src/clock/reroot.rs#L152-L184). Regression tests verify that the public reroot operation selects the same edge and split as a direct fixed-zero-rate search and that changing sampling-date values does not change a `min-dev` result [packages/treetime/src/clock/__tests__/test_reroot.rs#L77-L132](../../packages/treetime/src/clock/__tests__/test_reroot.rs#L77-L132).

This is an intentional correction of a v0 defect, not behavioral parity with v0's executed score. The `least-squares` mode preserves the estimated-rate temporal objective.

The approved v1 behavior is recorded in [kb/decisions/clock-min-dev-fixed-zero-rate.md](../decisions/clock-min-dev-fixed-zero-rate.md).

## Glossary

1. <a id="gloss-1"></a> **Root-to-tip distance.** The sum of branch lengths from a proposed tree root to a sampled leaf. [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Weighted least squares (WLS).** Regression that gives each observation a precision weight, so observations with lower modeled variance contribute more strongly. [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **Heterochronous samples.** Sequences collected at evolutionarily distinct times, allowing genetic divergence to be related to elapsed time. [↩](#gloss-use-3)
4. <a id="gloss-4"></a> **Strict molecular clock.** A model in which every branch accumulates substitutions at one common rate per unit time. [↩](#gloss-use-4)

## References

1. <a id="ref-1"></a> Mai, Uyen, Erfan Sayyari, and Siavash Mirarab. 2017. "Minimum Variance Rooting of Phylogenetic Trees and Implications for Species Tree Reconstruction." _PLOS ONE_ 12 (8): e0182238. https://doi.org/10.1371/journal.pone.0182238 [↩¹](#cite-1a) [↩²](#cite-1b)
2. <a id="ref-2"></a> Rambaut, Andrew, Tommy T. Lam, Luiz Max Carvalho, and Oliver G. Pybus. 2016. "Exploring the Temporal Structure of Heterochronous Sequences Using TempEst (Formerly Path-O-Gen)." _Virus Evolution_ 2 (1): vew007. https://doi.org/10.1093/ve/vew007 [↩](#cite-2)
3. <a id="ref-3"></a> Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4 (1): vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-3)
