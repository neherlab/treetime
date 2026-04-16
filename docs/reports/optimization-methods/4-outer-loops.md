# Chapter 4: Outer loops and EM convergence

[Back to index](_index.md) | Previous: [Chapter 3: Coalescent skyline](3-coalescent-skyline.md) | Next: [Chapter 5: Supporting optimizations](5-supporting-optimizations.md)

The per-edge branch length optimization (Chapter 2) runs inside an outer loop that alternates between ancestral reconstruction and parameter updates. This outer loop has an EM (Expectation-Maximization) structure. See [Iterative tree refinement, Chapter 9](../iterative-tree-refinement/9-iteration-loop.md) for the pipeline-level view.

## EM structure in phylogenetics

The phylogenetic likelihood has ancestral states at internal nodes as latent variables. The pruning algorithm marginalizes over them. When branch length optimization is cast as EM:

- **E-step**: compute posterior distribution of ancestral states (marginal reconstruction via backward + forward pruning pass)
- **M-step**: update branch lengths given expected sufficient statistics (per-edge optimization)

<a id="cite-1"></a>[Dempster, Laird, and Rubin 1977](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x) [[1](#ref-1)] established the EM algorithm with monotone likelihood increase guarantee. <a id="cite-2"></a>[Wu 1983](https://projecteuclid.org/journals/annals-of-statistics/volume-11/issue-1/On-the-Convergence-Properties-of-the-EM-Algorithm/10.1214/aos/1176346060.full) [[2](#ref-2)] corrected convergence proofs and extended beyond exponential families.

<a id="cite-3"></a>[Meng and Rubin 1993](https://doi.org/10.1093/biomet/80.2.267) [[3](#ref-3)] introduced ECM (Expectation-Conditional Maximization): replacing the single M-step with sequential conditional maximizations while preserving monotone convergence. TreeTime's alternating reconstruction + branch optimization + damping follows this ECM pattern.

## Convergence rate

EM convergence is linear (first-order), governed by the fraction of missing information - the ratio of information in the complete data to observed data (<a id="cite-4"></a>[Meng and Rubin 1994](<https://doi.org/10.1016/0024-3795(94)90363-8>) [[4](#ref-4)]). When sequences are long and the tree is shallow (TreeTime's viral phylogenetics use case), little information is lost by treating ancestral states as latent, and convergence is fast (2-5 iterations). For deep divergence or short sequences, convergence degrades.

<a id="cite-5"></a>[Sullivan et al. 2005](https://doi.org/10.1080/10635150590905966) [[5](#ref-5)] validated that alternating model parameter optimization and topology+branch-length optimization ("successive approximations") is reliable for ML topology estimation.

## Damping

TreeTime v1 uses exponential damping:

```
bl_damped = bl_new * (1 - d^(i+1)) + bl_old * d^(i+1)
```

with d=0.75. Early iterations are conservative (iter 0: 25% new, 75% old), later iterations are aggressive (iter 10: 94% new). This is a form of under-relaxation analogous to SOR from linear algebra. It preserves the monotone likelihood guarantee because the damped value is a convex combination of the old (known-good) and new (improved) branch lengths.

v1 code: [packages/treetime/src/commands/optimize/run.rs#L199-L221](../../../packages/treetime/src/commands/optimize/run.rs#L199-L221) `apply_damping()`

## How tools handle the outer loop

No tool uses pure EM for branch lengths. All use coordinate ascent (one branch at a time) with NR or Brent, iterated over all branches:

- **RAxML-NG**: up to 32 smoothings, adaptive epsilon
- **IQ-TREE**: up to 100 rounds, rollback if likelihood drops > tolerance\*0.1
- **PhyML**: until dLH < tolerance, alternating with model parameter optimization

## Acceleration methods

EM acceleration methods exist but are not adopted by any major phylogenetic tool:

- **SQUAREM** (<a id="cite-6"></a>[Varadhan and Roland 2008](https://doi.org/10.1111/j.1467-9469.2007.00585.x) [[6](#ref-6)]): off-the-shelf acceleration achieving superlinear convergence. Mean 18-fold speedup in benchmarks. Requires only the EM mapping function (no gradients).
- **DAAREM** ([Henderson and Varadhan 2019](https://doi.org/10.1080/10618600.2019.1594835)): damped Anderson acceleration with monotonicity control.
- **Anderson acceleration**: mixing of multiple previous iterates.

Whether these would accelerate TreeTime's outer loop is unknown. The primary bottleneck is typically the pruning traversal (O(n \* s^2) per site per pass), not the number of outer iterations.

## References

1. <a id="ref-1"></a> Dempster, Arthur P., Nan M. Laird, and Donald B. Rubin. 1977. "Maximum Likelihood from Incomplete Data Via the EM Algorithm." _JRSS:B_ 39(1):1-38. https://doi.org/10.1111/j.2517-6161.1977.tb01600.x [↩](#cite-1)
2. <a id="ref-2"></a> Wu, C. F. Jeff. 1983. "On the Convergence Properties of the EM Algorithm." _Ann. Stat._ 11(1):95-103. https://projecteuclid.org/journals/annals-of-statistics/volume-11/issue-1/On-the-Convergence-Properties-of-the-EM-Algorithm/10.1214/aos/1176346060.full (DOI: https://doi.org/10.1214/aos/1176346060) [↩](#cite-2)
3. <a id="ref-3"></a> Meng, Xiao-Li, and Donald B. Rubin. 1993. "Maximum Likelihood Estimation via the ECM Algorithm: A General Framework." _Biometrika_ 80(2):267-278. https://doi.org/10.1093/biomet/80.2.267 [↩](#cite-3)
4. <a id="ref-4"></a> Meng, Xiao-Li, and Donald B. Rubin. 1994. "On the Global and Componentwise Rates of Convergence of the EM Algorithm." _Lin. Alg. Appl._ 199:413-425. https://doi.org/10.1016/0024-3795(94)90363-8 [↩](#cite-4)
5. <a id="ref-5"></a> Sullivan, Jack, Zaid Abdo, Peter Joyce, and David L. Swofford. 2005. "Evaluating the Performance of a Successive-Approximations Approach to Parameter Optimization in Maximum-Likelihood Phylogeny Estimation." _Mol. Biol. Evol._ 22(6):1386-1392. https://doi.org/10.1080/10635150590905966 [↩](#cite-5)
6. <a id="ref-6"></a> Varadhan, Ravi, and Christophe Roland. 2008. "Simple and Globally Convergent Methods for Accelerating the Convergence of Any EM Algorithm." _Scand. J. Stat._ 35(2):335-353. https://doi.org/10.1111/j.1467-9469.2007.00585.x [↩](#cite-6)
