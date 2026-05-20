# Should indel rate be hoisted outside the optimize loop?

## The proposal

Hoist `estimate_indel_rate()` [packages/treetime/src/optimize/indel.rs#L55](../../packages/treetime/src/optimize/indel.rs#L55) out of the optimize loop, computing it once before iteration begins. The reasoning: the GTR substitution rate `gtr.mu` is estimated once and held fixed, so the indel rate should follow the same pattern.

## Why re-estimating is the right call

The analogy to `gtr.mu` is reasonable at first glance - both are global rate parameters. But the indel rate has a direct, explicit dependence on branch lengths that `gtr.mu` does not.

The indel rate MLE is $\hat{\mu} = \sum_e k_e \;/\; \sum_e t_e$. Branch lengths $t_e$ sit in the denominator and change every iteration. When the optimizer doubles a branch length, the indel rate estimate halves. This is not a subtle second-order effect; it is the dominant term.

`gtr.mu` is estimated from sequence composition and rate matrix structure. Branch lengths enter only indirectly through normalization. In practice, re-estimating `gtr.mu` after branch-length updates would produce near-identical values. Re-estimating the indel rate can produce meaningfully different values.

Re-estimating each iteration is coordinate ascent: optimize branch lengths holding $\mu$ fixed, then update $\mu$ at the new branch lengths, repeat. This is a standard approach for jointly dependent parameters.

### ECM justification

The optimize loop is an <a id="gloss-use-1"></a>ECM <sup>[1](#gloss-1)</sup> algorithm (<a id="cite-1"></a>[Meng and Rubin 1993](https://doi.org/10.1093/biomet/80.2.267) [[1](#ref-1)]). The latent data is the set of ancestral sequences. The <a id="gloss-use-2"></a>EM <sup>[2](#gloss-2)</sup> steps map as follows:

- E-step: marginal ancestral reconstruction (`update_marginal`) estimates ancestral sequences given current branch lengths
- CM-step 1: per-edge branch-length optimization (`run_optimize_mixed_with_indel_rate`) maximizes likelihood given the reconstructed ancestral sequences, holding the indel rate fixed
- CM-step 2: re-estimate the indel rate $\hat{\mu}$ at the new branch lengths (`estimate_indel_rate`)

ECM preserves EM's monotone convergence guarantee as long as every CM-step is executed each iteration. Freezing the indel rate permanently skips CM-step 2, turning the loop into a different (fixed-rate) optimization problem.

### Phylogenetic precedent

All major ML phylogenetic software re-estimates model parameters jointly with branch lengths in each optimization round. <a id="cite-2a"></a>[Guindon and Gascuel 2003](https://doi.org/10.1080/10635150390235520) [[2](#ref-2)] describe PhyML's hill-climbing algorithm as adjusting topology and branch lengths simultaneously, with model parameters updated per round. <a id="cite-3"></a>[Kozlov et al. 2019](https://doi.org/10.1093/bioinformatics/btz305) [[3](#ref-3)] note that RAxML-NG supports either optimizing or fixing all model parameters including branch lengths, with the default being joint optimization. <a id="cite-4"></a>[Nguyen et al. 2015](https://doi.org/10.1093/molbev/msu300) [[4](#ref-4)] describe IQ-TREE's stochastic perturbation method with per-round parameter updates. The standard practice since <a id="cite-5"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[5](#ref-5)] is to re-estimate all free parameters each round, not to freeze subsets.

The fact that `gtr.mu` is currently fixed in our codebase is a simplification, not a principled modeling choice. The convergence proposal [../proposals/optimize-convergence-and-robustness.md](../proposals/optimize-convergence-and-robustness.md) already identifies GTR re-estimation as a missing feature (P9), matching v0's `infer_gtr` parameter and the behavior of <a id="cite-2b"></a>[PhyML](https://doi.org/10.1080/10635150390235520) [[2](#ref-2)], RAxML-NG, and IQ-TREE.

## What would change if we hoisted it

The indel rate would reflect initial branch lengths only. As the optimizer rescales the tree, the Poisson indel penalty would use a stale rate, over- or under-penalizing indel-bearing edges relative to the current scale.

Whether this matters in practice depends on how much total branch length changes during optimization. If branch lengths barely move, the effect is negligible. If the tree rescales substantially (common in early iterations), the stale rate could bias convergence.

### Empirical evidence from the codebase

The intentional changes document [../decisions/optimize-indel-contribution-to-likelihood.md](../decisions/optimize-indel-contribution-to-likelihood.md) records that per-iteration $\hat\mu$ recomputation amplifies a 2-cycle on sc2/2844 ($\hat\mu \approx 12{,}000$, oscillating proportionally with total branch length). This is a real convergence concern, but the root cause is the sparse variable/fixed position boundary (not the indel rate itself). The convergence proposal's [P4](../proposals/optimize-convergence-and-robustness.md#L100) suggests caching $\hat\mu$ once topology stabilizes as a pragmatic fix, while P1/P2/P3 address the root cause.

Caching the indel rate is a valid workaround for the 2-cycle, but it changes the optimization objective. The clean fix is to address the sparse boundary discontinuity (P1-P3), which eliminates the 2-cycle without freezing a parameter that should be free.

## Summary

- Re-estimating each iteration is the theoretically correct ECM approach
- Standard phylogenetic ML software re-estimates all model parameters per round
- The `gtr.mu` precedent in this codebase is a known missing feature (P9), not a design principle to extend
- The 2-cycle amplification on sc2/2844 is real but caused by the sparse boundary, not by indel rate re-estimation
- If caching is adopted as a workaround, it should be documented as a model change with the intent to revert once P1-P3 resolve the root cause

## Glossary

1. <a id="gloss-1"></a> **ECM (Expectation Conditional Maximization).** EM variant that replaces the M-step with several conditional maximization (CM) steps, each optimizing one parameter subset while holding the rest fixed. Preserves EM's monotone convergence guarantee as long as every subset is updated each iteration (<a id="cite-1b"></a>[Meng and Rubin 1993](https://doi.org/10.1093/biomet/80.2.267) [[1](#ref-1)]). [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **EM (Expectation-Maximization).** Iterative algorithm for maximum-likelihood estimation when some data is unobserved (latent). E-step: estimate latent data using current parameters. M-step: re-estimate all parameters by maximizing likelihood given the estimated latent data. Guaranteed to increase likelihood each iteration (monotone convergence). [↩](#gloss-use-2)

## References

1. <a id="ref-1"></a> Meng, Xiao-Li, and Donald B. Rubin. 1993. "Maximum Likelihood Estimation via the ECM Algorithm: A General Framework." _Biometrika_ 80:267-278. https://doi.org/10.1093/biomet/80.2.267 [↩¹](#cite-1) [↩²](#cite-1b)
2. <a id="ref-2"></a> Guindon, Stéphane, and Olivier Gascuel. 2003. "A Simple, Fast, and Accurate Algorithm to Estimate Large Phylogenies by Maximum Likelihood." _Systematic Biology_ 52:696-704. https://doi.org/10.1080/10635150390235520 [↩¹](#cite-2a) [↩²](#cite-2b)
3. <a id="ref-3"></a> Kozlov, Alexey M., Diego Darriba, Tomáš Flouri, Benoit Morel, and Alexandros Stamatakis. 2019. "RAxML-NG: A Fast, Scalable and User-Friendly Tool for Maximum Likelihood Phylogenetic Inference." _Bioinformatics_ 35:4453-4455. https://doi.org/10.1093/bioinformatics/btz305 [↩](#cite-3)
4. <a id="ref-4"></a> Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh. 2015. "IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies." _Molecular Biology and Evolution_ 32:268-274. https://doi.org/10.1093/molbev/msu300 [↩](#cite-4)
5. <a id="ref-5"></a> Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17:368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-5)
