# Per-iteration indel rate re-estimation

## Proposal

Revert the hoisting of `estimate_indel_rate()` from before the optimize loop to inside each iteration, restoring per-iteration re-estimation as a coordinate-ascent CM-step.

## Current state

`estimate_indel_rate()` is computed once before the first iteration in `run_optimize_loop()` and held fixed throughout. This removes the feedback path where branch-length changes shift $\hat\mu = \sum_e k_e / \sum_e t_e$, which shifts branch-length targets. The `--no-indels` flag separately disables the indel contribution entirely.

## Motivation

Per-iteration re-estimation is the theoretically correct ECM (Expectation Conditional Maximization) approach. The optimize loop is an ECM algorithm where:

- E-step: marginal ancestral reconstruction estimates ancestral sequences given current branch lengths
- CM-step 1: per-edge branch-length optimization maximizes likelihood given reconstructed ancestral sequences, holding the indel rate fixed
- CM-step 2: re-estimate the indel rate $\hat\mu$ at the new branch lengths

ECM preserves EM's monotone convergence guarantee as long as every CM-step executes each iteration ([Meng and Rubin 1993](https://doi.org/10.1093/biomet/80.2.267)). Permanently freezing the indel rate skips CM-step 2, turning the loop into a different (fixed-rate) optimization problem.

### Phylogenetic precedent

All major ML phylogenetic software re-estimates model parameters jointly with branch lengths per optimization round: [PhyML](https://doi.org/10.1080/10635150390235520) (Guindon and Gascuel 2003), [RAxML-NG](https://doi.org/10.1093/bioinformatics/btz305) (Kozlov et al. 2019), [IQ-TREE](https://doi.org/10.1093/molbev/msu300) (Nguyen et al. 2015). The standard practice since [Felsenstein 1981](https://doi.org/10.1007/BF01734359) is to re-estimate all free parameters each round, not to freeze subsets.

### The `gtr.mu` precedent is not a precedent

The analogy to `gtr.mu` (estimated once and held fixed) is reasonable at first glance but misleading. The indel rate has a direct, explicit dependence on branch lengths ($t_e$ sits in the denominator of $\hat\mu$), while `gtr.mu` is estimated from sequence composition and rate matrix structure with only indirect branch-length dependence. Re-estimating `gtr.mu` after branch-length updates would produce near-identical values; re-estimating the indel rate can produce meaningfully different values.

The fact that `gtr.mu` is currently fixed is itself a known missing feature documented as P9 in the [convergence and robustness proposal](optimize-convergence-and-robustness.md), not a design principle to extend.

### 2-cycle amplification root cause

On sc2/2844, $\hat\mu \approx 12{,}000$ oscillates proportionally with total branch length across iterations. Per-iteration re-estimation amplifies this 2-cycle, but the root cause is the sparse variable/fixed position boundary discontinuity (P1-P3 in the convergence proposal), not the indel rate re-estimation itself. The clean fix is to address the boundary discontinuity, which eliminates the 2-cycle without freezing a parameter that should be free.

## Implementation

Restore `estimate_indel_rate()` inside the loop body of `run_optimize_loop()`, computing a fresh $\hat\mu$ each iteration after `update_marginal()` and before `run_optimize_mixed_inner()`. The `--no-indels` gate remains: when set, the rate stays at 0.0 regardless.

## Validation

- Existing golden master tests verify branch-length parity with v0 (v0 does not have an indel contribution, so this change is orthogonal to v0 parity)
- Monitor convergence behavior on sc2/2844 for 2-cycle amplification. If the 2-cycle returns, P1-P3 from the convergence proposal should be prioritized rather than re-freezing the rate

## Dependencies

- P1-P3 from [convergence and robustness proposal](optimize-convergence-and-robustness.md) address the root cause of the 2-cycle that motivated the hoisting
- P9 from the same proposal addresses `gtr.mu` re-estimation, aligning it with the indel rate treatment

## Full analysis

See [indel rate re-estimation report](../reports/optimize-indel-rate-reestimation.md) for the complete ECM justification, phylogenetic precedent, and empirical evidence.

## References

- Meng, Xiao-Li, and Donald B. Rubin. 1993. "Maximum Likelihood Estimation via the ECM Algorithm: A General Framework." _Biometrika_ 80(2):267-278. https://doi.org/10.1093/biomet/80.2.267
- Guindon, Stephane, and Olivier Gascuel. 2003. "A Simple, Fast, and Accurate Algorithm to Estimate Large Phylogenies by Maximum Likelihood." _Systematic Biology_ 52(5):696-704. https://doi.org/10.1080/10635150390235520
- Kozlov, Alexey M., Diego Darriba, Tomas Flouri, Benoit Morel, and Alexandros Stamatakis. 2019. "RAxML-NG: A Fast, Scalable and User-Friendly Tool for Maximum Likelihood Phylogenetic Inference." _Bioinformatics_ 35(21):4453-4455. https://doi.org/10.1093/bioinformatics/btz305
- Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh. 2015. "IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies." _Molecular Biology and Evolution_ 32(1):268-274. https://doi.org/10.1093/molbev/msu300
- Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359
