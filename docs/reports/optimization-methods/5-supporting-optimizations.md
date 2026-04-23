# Chapter 5: Supporting optimizations

[Back to index](_index.md) | Previous: [Chapter 4: Outer loops](4-outer-loops.md) | Next: [Chapter 6: The argmin crate](6-argmin-crate.md)

Beyond branch lengths (Chapter 2) and skyline (Chapter 3), TreeTime performs several supporting optimizations. Each has different structure and optimizer requirements.

## Root-split optimization

Finding the optimal position along an edge to place the tree root. Minimizes chi-squared of the clock regression as a function of position `x in [0, 1]`. No analytical derivatives are available.

v1 provides three interchangeable solvers (E1-E3): BrentOpt, GoldenSectionSearch, and grid search. All use the same `BranchPointCostFunction`. See [audit inconsistency I1](7-audit.md) for why grid search should be retired.

v0 uses `scipy.optimize.minimize_scalar(method='bounded')` with a 6-point grid pre-search: [treetime/treeregression.py#L385-L414](https://github.com/neherlab/treetime/blob/master/treetime/treeregression.py#L385-L414).

## Coalescent Tc optimization

Finding the constant coalescent time scale Tc that maximizes the total coalescent likelihood. BrentOpt in log-space with bracket [-20, 2]. Already uses argmin (E4).

## Polytomy resolution

Greedy pairwise merging of children under new internal nodes. For each candidate pair, BrentOpt finds the optimal merge time (E6). The outer loop picks the best pair and repeats until no pair exceeds the resolution threshold (default 0.05).

v0 deprecated greedy resolution in favor of stochastic coalescent resolution (`generate_subtree()`), which produces coalescent-typical tree shapes instead of caterpillar-biased ones. v1 has not ported the stochastic method (`N-timetree-stochastic-polytomy-unimplemented`).

<a id="cite-1"></a>[Lewis, Holder, and Holsinger 2005](https://doi.org/10.1080/10635150590924208) [[1](#ref-1)] showed that Bayesian MCMC assigns near-1.0 posterior probability to arbitrary binary resolutions of true star trees (the star tree paradox). TreeTime's `resolution_threshold` serves an analogous role to their polytomy prior: preventing commitment to unsupported resolutions.

## HPD interval computation

Finding the highest posterior density region containing a specified fraction of probability mass. The algorithm bisects on a density threshold `p_thresh in [0, peak_val]` to find where the enclosed mass equals the target fraction.

This is a 1D root-finding problem. The current hand-rolled bisection converges in ~50 iterations to machine precision. argmin's `BrentRoot` would converge in ~10-15 iterations but the function is called once per node with O(1) cost per iteration, so the benefit is marginal.

<a id="cite-2"></a>[Hyndman 1996](https://doi.org/10.1080/00031305.1996.10474359) [[2](#ref-2)] established the threshold-based framework. <a id="cite-3"></a>[Chen and Shao 1999](https://doi.org/10.1080/10618600.1999.10474802) [[3](#ref-3)] developed the sorted-sample method used in BEAST/Tracer.

## GTR parameter estimation

Inferring substitution model parameters (rate matrix W, frequencies pi, rate mu) from observed mutation counts. TreeTime uses a fixed-point iteration: alternating closed-form updates of W, pi, and mu until convergence. This is an ECM scheme similar to <a id="cite-4"></a>[Holmes and Rubin 2002](https://doi.org/10.1006/jmbi.2002.5405) [[4](#ref-4)] but uses MAP ancestral reconstruction as a hard assignment rather than full posterior expectations.

Other tools use different approaches:

- RAxML-NG, PAML: quasi-Newton (BFGS) joint optimization of all rate + frequency parameters
- PhyML: cyclic coordinate descent (Brent per parameter)
- XRate: EM with full posterior sufficient statistics

argmin is not suitable for TreeTime's fixed-point iteration (audit item H2).

## Relaxed clock

Computing per-branch rate multipliers (gamma) via an autocorrelated model with quadratic penalty. The penalty has three terms: slack (deviation from mean rate), stiffness (mismatch between rescaled branch length and observed mutations), and coupling (rate difference between parent and child).

The two-pass algorithm (postorder coefficient accumulation, preorder gamma computation) solves the optimization in O(n) time with no iteration. Each gamma is computed analytically as `max(0.1, -0.5 * k1 / k2)`.

TreeTime's coupling penalty `(gamma_parent - gamma_child)^2` is the Gaussian analog of the lognormal prior in <a id="cite-5"></a>[Thorne, Kishino, and Painter 1998](https://doi.org/10.1093/oxfordjournals.molbev.a025892) [[5](#ref-5)]. Setting coupling=0 recovers an uncorrelated model (<a id="cite-6"></a>[Drummond et al. 2006](https://doi.org/10.1371/journal.pbio.0040088) [[6](#ref-6)]).

## References

1. <a id="ref-1"></a> Lewis, Paul O., Mark T. Holder, and Kent E. Holsinger. 2005. "Polytomies and Bayesian Phylogenetic Inference." _Syst. Biol._ 54(2):241-253. https://doi.org/10.1080/10635150590924208 [↩](#cite-1)
2. <a id="ref-2"></a> Hyndman, Rob J. 1996. "Computing and Graphing Highest Density Regions." _The American Statistician_ 50(2):120-126. https://doi.org/10.1080/00031305.1996.10474359 [↩](#cite-2)
3. <a id="ref-3"></a> Chen, Ming-Hui, and Qi-Man Shao. 1999. "Monte Carlo Estimation of Bayesian Credible and HPD Intervals." _J. Comput. Graph. Stat._ 8(1):69-92. https://doi.org/10.1080/10618600.1999.10474802 [↩](#cite-3)
4. <a id="ref-4"></a> Holmes, Ian, and Gregory M. Rubin. 2002. "An Expectation Maximization Algorithm for Training Hidden Substitution Models." _J. Mol. Biol._ 317(5):753-764. https://doi.org/10.1006/jmbi.2002.5405 [↩](#cite-4)
5. <a id="ref-5"></a> Thorne, Jeffrey L., Hirohisa Kishino, and Ian S. Painter. 1998. "Estimating the Rate of Evolution of the Rate of Molecular Evolution." _Mol. Biol. Evol._ 15(12):1647-1657. https://doi.org/10.1093/oxfordjournals.molbev.a025892 [↩](#cite-5)
6. <a id="ref-6"></a> Drummond, Alexei J., Simon Y. W. Ho, Matthew J. Phillips, and Andrew Rambaut. 2006. "Relaxed Phylogenetics and Dating with Confidence." _PLoS Biology_ 4(5):e88. https://doi.org/10.1371/journal.pbio.0040088 [↩](#cite-6)
