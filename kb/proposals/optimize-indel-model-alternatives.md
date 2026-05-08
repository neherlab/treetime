# Proposal: indel model alternatives for branch length optimization

## Summary

Evaluate whether the current Poisson indel count model should be extended with length-aware weighting or replaced with a richer model for datasets where indel events carry substantial phylogenetic signal.

## Current state

v1 uses a Poisson count model where each indel event contributes equally to the branch length log-likelihood regardless of indel length. This prevents zero-length assignment on branches with only indel evidence. See [intentional change](../decisions/optimize-indel-contribution-to-likelihood.md) and [indel models report](../reports/indel-models/README.md) for the full catalog of indel models and the rationale for the current choice.

The model is implemented in [`packages/treetime/src/commands/optimize/optimize_indel.rs`](../../packages/treetime/src/commands/optimize/optimize_indel.rs).

## Limitation

The equal-weight assumption means a 100-base deletion contributes the same likelihood as a 1-base deletion. For typical viral phylodynamics (TreeTime's primary use case), indels are rare and short, so this has negligible impact. For datasets with frequent, long indels (bacterial genomics, gene family evolution, recombination-prone viruses), the model underweights long indel events.

## Candidate extensions

The extensions below are ordered by implementation complexity. Each preserves the additive structure with the substitution log-likelihood in the Newton step. None require restructuring the per-edge likelihood evaluation.

### E1. Length-weighted Poisson count

Replace the indel count $k$ with a length-weighted count $k_w = \sum_i w(l_i)$ where $l_i$ is the length of the $i$-th indel and $w$ is a weighting function. The Poisson log-likelihood becomes:

$$\ell_{\text{indel}}(t) = k_w \ln(\mu_w t) - \mu_w t - \text{const}$$

where $\mu_w = \sum_e k_{w,e} / \sum_e t_e$.

Weighting options:

- Linear: $w(l) = l$. A 100-base deletion contributes 100x as much as a 1-base deletion. The indel rate $\mu_w$ becomes "total indel bases per unit branch length."
- Logarithmic: $w(l) = \ln(l + 1)$. Diminishing returns for longer indels. Motivated by the observation that long indels are rarer but individually more informative.
- Square root: $w(l) = \sqrt{l}$. Intermediate between count and linear.

The factorial term $\ln(k!)$ in the exact Poisson is not well-defined for non-integer $k_w$. Use the Stirling approximation or drop it (it is $t$-independent and does not affect the Newton step).

**Implementation effort**: minimal. The `edge_indel_count()` trait method already returns `usize`; extend to return a weighted sum. The `InDel` struct stores position and length, so length information is available.

**Tradeoff**: adds no new parameters. The choice of weighting function is a modeling decision, not a fitted parameter. Linear weighting is the simplest and most interpretable.

### E2. Separate insertion and deletion rates

Replace the single rate $\mu$ with separate insertion rate $\lambda$ and deletion rate $\delta$:

$$\ell_{\text{indel}}(t) = k_{\text{ins}} \ln(\lambda t) - \lambda t + k_{\text{del}} \ln(\delta t) - \delta t - \text{const}$$

<a id="cite-1"></a>[Loewenthal et al. 2021](https://doi.org/10.1093/molbev/msab266) [[1](#ref-1)] found 74% of RIM-supported protein datasets have $R_D > R_I$, with $R_D/R_I \approx 2$ in Drosophilidae. The asymmetry is not universal: Saccharomycetaceae coding genes show $R_I > R_D$ in 56% of RIM datasets.

**Mathematical note**: with global MLE estimation, $\lambda + \delta = \mu$ always holds (the total rate is the same). The per-edge optimum $t^* = (k_{\text{ins}} + k_{\text{del}})/(\lambda + \delta) = k/\mu$ is identical to the current single-rate model. Separate rates do not change branch length estimates under global MLE. The separation is useful for reporting and for downstream analyses that distinguish insertion-dominated from deletion-dominated branches, but it requires per-edge or per-clade rate estimation to affect branch lengths.

**Tradeoff**: no branch length change with global rates. Requires per-edge or hierarchical rate estimation to produce different branch lengths, which adds complexity.

### E3. Affine-inspired two-parameter model

Model the indel contribution as a sum of an opening cost and an extension cost, analogous to affine gap penalties but within a probabilistic framework:

$$\ell_{\text{indel}}(t) = k \ln(\mu_o t) + L_{\text{ext}} \ln(\mu_e t) - (\mu_o + \mu_e) t$$

where $k$ is the number of indel events, $L_{\text{ext}} = \sum_i (l_i - 1)$ is the total extension length, $\mu_o$ is the gap opening rate, and $\mu_e$ is the gap extension rate.

This factorizes the indel process into "how many events" (opening) and "how long each event" (extension), matching the biological intuition that gap initiation and gap propagation have different evolutionary dynamics.

**Implementation effort**: moderate. Two rate parameters estimated from the tree. The Newton step gains two terms. Requires the `InDel` struct to store length (already available).

**Mathematical note**: the per-edge optimum $t^* = (k + L_{\text{ext}})/(\mu_o + \mu_e)$. Since $k + L_{\text{ext}} = \sum l_i$, E3 is equivalent to E1 with linear weighting $w(l) = l$ when $\mu_o = \mu_e$. The separate opening/extension rates add a degree of freedom beyond E1: they allow the model to weight short indels (dominated by opening) differently from long indels (dominated by extension). <a id="cite-2"></a>[Rivas 2005](https://doi.org/10.1186/1471-2105-6-236) [[2](#ref-2)] showed geometric instantaneous rates do not produce geometric finite-time distributions, so this is an approximation.

**Tradeoff**: two parameters. More expressive than E1 (separates opening from extension dynamics). Subsumes E1-linear as a special case ($\mu_o = \mu_e$).

## Architectural alternatives

These require deeper changes to the optimization architecture. Each is described with tradeoffs for evaluation.

### E4. TKF91/TKF92 integration

Integrate TKF91 or TKF92 as the indel model for branch length optimization. TKF91 jointly models substitutions and indels via a pair HMM with separate insertion rate $\lambda$ and deletion rate $\mu$. Branch length optimization uses Brent's method (as in rust-phylo) or Newton-Raphson with pair HMM derivatives.

**Fixed-alignment mode**: with a fixed alignment (treating it as observed data), TKF91 likelihood per edge is $O(L)$ via column-by-column evaluation. This is the mode relevant to TreeTime's optimizer, where the alignment is fixed and only branch lengths are optimized. The $O(L^N)$ cost applies to marginalizing over all possible alignments, which TreeTime does not do.

**Architectural change**: replaces the per-site eigendecomposition coefficient cache with a per-column pair HMM evaluation. The rust-phylo library (`acg-team/rust-phylo`) demonstrates this architecture in Rust with `argmin::BrentOpt` for per-edge optimization.

**Tradeoff**: 2 parameters ($\lambda$, $\mu$). Distinguishes insertions from deletions. Single-residue indels only (TKF92 extends to geometric fragments with parameter $r$). Requires pair HMM DP per edge instead of scalar addition.

### E5. Poisson Indel Process (PIP)

PIP decouples insertions from existing characters, making the marginal likelihood (summing over alignments) linear in the number of taxa. For fixed-alignment branch length optimization, PIP provides a column-based likelihood where each column's contribution depends on whether the character was inserted or deleted on each branch.

**Architectural change**: requires tracking per-column insertion/deletion status on each edge, not just a total count. The column-based likelihood replaces the scalar Poisson term with a per-column product.

**Tradeoff**: 2 parameters ($\lambda$, $\mu$). Linear-time marginal likelihood (relevant if TreeTime ever marginalizes over alignments). Equilibrium length distribution is Poisson (lighter tails than TKF91's geometric).

### E6. Length-distribution fitting

Fit indel length distributions (geometric, Zipf/power-law, multi-exponential) from the observed indel data, then use the fitted distribution to weight each indel event in the branch length likelihood.

**Approach**: estimate distribution parameters from the tree's indel length histogram (e.g., MLE for geometric, or ABC as in SpartaABC for Zipf). Use the fitted probability $P(l_i)$ as a weight for each indel event in the Poisson term: $\ell(t) = \sum_i \ln(P(l_i)) + k \ln(\mu t) - \mu t$.

**Tradeoff**: 1-2 additional distribution parameters. Requires a nested optimization (fit distribution, then optimize branch lengths, iterate). More principled than E1's fixed weighting function but adds iteration complexity.

## Validation plan

1. Compare branch lengths with and without the extension on datasets with varying indel frequencies:
   - Low indel: flu/h3n2/20, ebola (viral, rare indels)
   - Moderate indel: tb, lassa/L (bacterial/viral with structural variation)
   - Synthetic: generated trees with known indel rates and lengths

2. For each extension, verify:
   - Log-likelihood at the optimum is higher than the current model
   - The optimizer converges (no oscillation, no NaN)
   - Branch lengths on indel-only edges are proportional to indel content
   - No regression on datasets where indels are absent (the extension should be a no-op)

3. Compare against v0 (which ignores indels entirely) to quantify the practical impact.

## References

- <a id="ref-1"></a>Loewenthal, Gil, Dana Rapoport, Oren Avram, Alon Moshe, Anat Kreimer, Michal Linial, and Tal Pupko. 2021. "A Probabilistic Model for Indel Evolution: Differentiating Insertions from Deletions." _Molecular Biology and Evolution_ 38(12):5769-5781. https://doi.org/10.1093/molbev/msab266 [↩](#cite-1)
- <a id="ref-2"></a>Rivas, Elena. 2005. "Evolutionary Models for Insertions and Deletions in a Probabilistic Modeling Framework." _BMC Bioinformatics_ 6:236. https://doi.org/10.1186/1471-2105-6-236 [↩](#cite-2)

## Related documents

- [Indel models report](../reports/indel-models/README.md) - full catalog of indel modeling approaches
- [Indel contribution intentional change](../decisions/optimize-indel-contribution-to-likelihood.md) - current Poisson model documentation
- ~~Grid search ignores indels~~ - **FIXED**
- ~~Timetree ignores indels~~ - **FIXED**
- [Convergence and robustness proposal](optimize-convergence-and-robustness.md) - related optimization proposal
