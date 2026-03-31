# Proposal: indel model alternatives for branch length optimization

## Summary

Evaluate whether the current Poisson indel count model should be extended with length-aware weighting or replaced with a richer model for datasets where indel events carry substantial phylogenetic signal.

## Current state

v1 uses a Poisson count model where each indel event contributes equally to the branch length log-likelihood regardless of indel length. This prevents zero-length assignment on branches with only indel evidence. See [intentional change](../port-intentional-changes/optimize-indel-contribution-to-likelihood.md) and [algorithm inventory](../port-algo-inventory/indel-models.md) for the full catalog of indel models and the rationale for the current choice.

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

- **Linear**: $w(l) = l$. A 100-base deletion contributes 100x as much as a 1-base deletion. The indel rate $\mu_w$ becomes "total indel bases per unit branch length."
- **Logarithmic**: $w(l) = \ln(l + 1)$. Diminishing returns for longer indels. Motivated by the observation that long indels are rarer but individually more informative.
- **Square root**: $w(l) = \sqrt{l}$. Intermediate between count and linear.

The factorial term $\ln(k!)$ in the exact Poisson is not well-defined for non-integer $k_w$. Use the Stirling approximation or drop it (it is $t$-independent and does not affect the Newton step).

**Implementation effort**: minimal. The `edge_indel_count()` trait method already returns `usize`; extend to return a weighted sum. The `InDel` struct stores position and length, so length information is available.

**Tradeoff**: adds no new parameters. The choice of weighting function is a modeling decision, not a fitted parameter. Linear weighting is the simplest and most interpretable.

### E2. Separate insertion and deletion rates

Replace the single rate $\mu$ with separate insertion rate $\lambda$ and deletion rate $\delta$:

$$\ell_{\text{indel}}(t) = k_{\text{ins}} \ln(\lambda t) - \lambda t + k_{\text{del}} \ln(\delta t) - \delta t - \text{const}$$

Empirical studies consistently find $\delta > \lambda$ (deletion rate exceeds insertion rate) across most datasets (Loewenthal et al. 2021, doi:10.1093/molbev/msab266).

**Implementation effort**: moderate. Requires distinguishing insertions from deletions in the `InDel` struct (already has `InDelType` enum with `Insert` and `Delete` variants). Rate estimation splits into two MLEs. The Newton step gains two additive terms instead of one.

**Tradeoff**: one additional parameter. Improves biological realism for datasets with asymmetric indel rates. For viral datasets where indels are rare, both rates are near zero and the distinction is immaterial.

### E3. Affine-inspired two-parameter model

Model the indel contribution as a sum of an opening cost and an extension cost, analogous to affine gap penalties but within a probabilistic framework:

$$\ell_{\text{indel}}(t) = k \ln(\mu_o t) + L_{\text{ext}} \ln(\mu_e t) - (\mu_o + \mu_e) t$$

where $k$ is the number of indel events, $L_{\text{ext}} = \sum_i (l_i - 1)$ is the total extension length, $\mu_o$ is the gap opening rate, and $\mu_e$ is the gap extension rate.

This factorizes the indel process into "how many events" (opening) and "how long each event" (extension), matching the biological intuition that gap initiation and gap propagation have different evolutionary dynamics.

**Implementation effort**: moderate. Two rate parameters estimated from the tree. The Newton step gains two terms. Requires the `InDel` struct to store length (already available).

**Tradeoff**: two parameters. The opening and extension rates are independently estimated, which is more principled than the affine gap penalty's fixed ratio. The model is still not a proper evolutionary process (Rivas 2005 showed geometric instantaneous rates do not produce geometric finite-time distributions), but it is a better approximation than equal-weight counting for branch length optimization.

## Not proposed

### TKF91/TKF92 integration

Full TKF91 or TKF92 integration would require restructuring the per-edge likelihood evaluation to use pair HMMs or extended Felsenstein peeling. The computational cost ($O(L^N)$ exact for TKF91) is prohibitive for TreeTime's interactive use case. The implementation effort is disproportionate to the benefit for branch length optimization specifically.

### Poisson Indel Process (PIP)

PIP achieves linear-time marginal likelihood by decoupling insertions from existing characters. While tractable, it still requires alignment-aware likelihood computation, not just a scalar indel count per edge. The integration would change the optimization architecture, not just add a term.

### Length-distribution fitting

Fitting indel length distributions (geometric, Zipf/power-law, multi-exponential) from the data would add statistical rigor but requires solving a nested optimization problem (estimate length distribution parameters, then use them in branch length optimization). The benefit for branch length optimization is marginal compared to simpler weighting schemes.

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

## Related documents

- [Indel models algorithm inventory](../port-algo-inventory/indel-models.md) - full catalog of indel modeling approaches
- [Indel contribution intentional change](../port-intentional-changes/optimize-indel-contribution-to-likelihood.md) - current Poisson model documentation
- [Design doc](../algorithms/optimize.md#poisson-indel-contribution-implemented) - model formulation and assumptions
- [Grid search ignores indels](../port-known-issues/M-optimize-grid-zero-ignores-indels.md) - related known issue
- [Timetree ignores indels](../port-known-issues/N-timetree-branch-length-distribution-ignores-indels.md) - related known issue
- [Convergence and method choice proposal](optimize-convergence-and-method-choice.md) - related optimization proposal
