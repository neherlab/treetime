# Clock min-dev uses a fixed zero clock rate

Clock `min-dev` ranks candidate roots with weighted least squares at a fixed clock rate of zero. It does not estimate a temporal slope from sampling dates.

**Type**: Bug fix (v0 erratum correction).

**v0 erratum**: [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md).

**v1 location**: `fn select_root()` dispatches `RerootMethod::MinDev` to `RootObjective::FixedRate(0.0)` [packages/treetime/src/clock/reroot.rs#L152-L184](../../packages/treetime/src/clock/reroot.rs#L152-L184).

## Decision

For candidate root $r$, clock `min-dev` minimizes

$$
J_0(r) = \min_{\alpha}\frac{1}{2W}\sum_{i=1}^{n}w_i\left(d_i(r)-\alpha\right)^2
= \frac{S_{dd}(r)W-S_d(r)^2}{2W^2}.
$$

Here, $J_0(r)$ is the fixed-zero-rate objective, $\alpha$ is the fitted intercept, $W=\sum_i w_i$ is total precision, $n$ is the number of contributing tips, $w_i$ is tip $i$'s precision, $d_i(r)$ is its distance from candidate root $r$, $S_d(r)=\sum_i w_i d_i(r)$, and $S_{dd}(r)=\sum_i w_i d_i(r)^2$.

This objective minimizes weighted root-to-tip variance, the minimum-variance rooting criterion described by <a id="cite-1"></a>[Mai et al. 2017](https://doi.org/10.1371/journal.pone.0182238) [[1](#ref-1)]. Sampling-date values are absent, so changing dates cannot alter the selected root when the contributing tips and variance model remain unchanged.

Clock `least-squares` remains the temporal-rooting mode. It estimates the clock rate and minimizes residuals around a root-to-tip regression, consistent with the methods described by <a id="cite-2"></a>[Rambaut et al. 2016](https://doi.org/10.1093/ve/vew007) [[2](#ref-2)] and <a id="cite-3"></a>[Sagulenko et al. 2018](https://doi.org/10.1093/ve/vex042) [[3](#ref-3)]. The two modes therefore answer distinct questions: uniformity of tip depths for `min-dev`, and temporal clock fit for `least-squares`.

## v0 divergence

TreeTime v0 documents `min_dev` as fixed-zero-rate minimum-variance rooting and supplies `slope=0`, but its candidate score still subtracts the residual reduction from an estimated temporal slope. v1 follows the documented and scientifically distinct objective instead of reproducing that scoring defect. The source-level evidence and full derivation are recorded in [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md).

## Validation contract

- The public clock reroot operation must select the same edge and split as a direct `RootObjective::FixedRate(0.0)` search.
- Changing only finite sampling-date values must leave the `min-dev` result unchanged.
- `least-squares` must continue to use `RootObjective::EstimatedRate`.

## References

1. <a id="ref-1"></a> Mai, Uyen, Erfan Sayyari, and Siavash Mirarab. 2017. "Minimum Variance Rooting of Phylogenetic Trees and Implications for Species Tree Reconstruction." _PLOS ONE_ 12 (8): e0182238. https://doi.org/10.1371/journal.pone.0182238 [↩](#cite-1)
2. <a id="ref-2"></a> Rambaut, Andrew, Tommy T. Lam, Luiz Max Carvalho, and Oliver G. Pybus. 2016. "Exploring the Temporal Structure of Heterochronous Sequences Using TempEst (Formerly Path-O-Gen)." _Virus Evolution_ 2 (1): vew007. https://doi.org/10.1093/ve/vew007 [↩](#cite-2)
3. <a id="ref-3"></a> Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4 (1): vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-3)
