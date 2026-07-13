# Algorithm Inventory - TreeTime v1 (Rust) and v0 (Python)

> **Valid for commit**: `688776a9` (rust branch, 2026-05-20)

Scope: algorithms implemented in v1 (domain files) and v0 algorithms not yet implemented ([unimplemented.md](unimplemented.md)). New v1 features not in v0 belong in [`proposals/`](../proposals/README.md).

## Summary

| Domain                                      | v1      | Unimplemented | Well-Known | Custom  |
| ------------------------------------------- | ------- | ------------- | ---------- | ------- |
| [Ancestral Reconstruction](ancestral.md)    | 2       | 1             | 3          | 0       |
| [Clock Inference](clock.md)                 | 9       | 3             | 8          | 8       |
| [Rerooting](reroot.md)                      | 2       | 0             | 0          | 2       |
| [Timetree Inference](timetree.md)           | 20      | 1             | 11         | 6       |
| [Mugration](mugration.md)                   | 3       | 1             | 2          | 0       |
| [Distribution/Convolution](distribution.md) | 17      | 4             | 12         | 5       |
| [GTR Substitution Models](gtr.md)           | 12      | 3             | 8          | 3       |
| [Graph Traversal](graph.md)                 | 4       | 0             | 4          | 0       |
| [Numerical Optimization](optimization.md)   | 14      | 1             | 5          | 0       |
| [Indel Models](indel-models.md)             | 1       | 8             | 7          | 1       |
| **Total**                                   | **~82** | **22**        | **~60**    | **~23** |

Full unimplemented algorithm details: [unimplemented.md](unimplemented.md)

---

## Quick Reference

| Algorithm                                                                       | Type       | Domain       | Location                                                                                                                                                                                            | Status        |
| ------------------------------------------------------------------------------- | ---------- | ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------- |
| [Fitch Parsimony](ancestral.md#fitch-parsimony)                                 | well-known | ancestral    | [`packages/treetime/src/ancestral/fitch.rs`](../../packages/treetime/src/ancestral/fitch.rs)                                                                                                        | complete      |
| [Marginal ML](ancestral.md#marginal-ml)                                         | well-known | ancestral    | [`packages/treetime/src/partition/marginal_*.rs`](../../packages/treetime/src/partition/)                                                                                                           | complete      |
| [Discrete Marginal](mugration.md#discrete-marginal-reconstruction)              | well-known | mugration    | [`packages/treetime/src/partition/marginal_core.rs`](../../packages/treetime/src/partition/marginal_core.rs) + [`marginal_discrete.rs`](../../packages/treetime/src/partition/marginal_discrete.rs) | complete      |
| [Joint ML](unimplemented.md#joint-ml)                                           | well-known | ancestral    | -                                                                                                                                                                                                   | unimplemented |
| [WLS Sufficient Stats](clock.md#wls-sufficient-statistics)                      | custom     | clock        | [`packages/treetime/src/payload/clock_set.rs`](../../packages/treetime/src/payload/clock_set.rs)                                                                                                    | complete      |
| [Tree Regression](clock.md#tree-regression)                                     | well-known | clock        | [`packages/treetime/src/clock/clock_regression.rs`](../../packages/treetime/src/clock/clock_regression.rs)                                                                                          | complete      |
| [Brent's Method](clock.md#brents-method)                                        | well-known | clock        | [`packages/treetime/src/clock/find_best_root/method_brent.rs`](../../packages/treetime/src/clock/find_best_root/method_brent.rs)                                                                    | complete      |
| [IQD Outlier Detection](clock.md#iqd-outlier-detection)                         | well-known | clock        | [`packages/treetime/src/clock/clock_filter.rs`](../../packages/treetime/src/clock/clock_filter.rs)                                                                                                  | complete      |
| [Generic Root Search](reroot.md#generic-root-search)                            | custom     | reroot       | [`packages/treetime/src/reroot/`](../../packages/treetime/src/reroot/)                                                                                                                              | complete      |
| [Divergence Rooting](reroot.md#divergence-rooting-min-dev)                      | custom     | reroot       | [`packages/treetime/src/reroot/div_stats.rs`](../../packages/treetime/src/reroot/div_stats.rs)                                                                                                      | complete      |
| [Belief Propagation](timetree.md#belief-propagation)                            | well-known | timetree     | [`packages/treetime/src/timetree/inference/`](../../packages/treetime/src/timetree/inference/)                                                                                                      | complete      |
| [Kingman Coalescent](timetree.md#kingman-coalescent)                            | well-known | timetree     | [`packages/treetime/src/coalescent/`](../../packages/treetime/src/coalescent/)                                                                                                                      | complete      |
| [Skyline Coalescent](timetree.md#skyline-coalescent)                            | well-known | timetree     | [`packages/treetime/src/coalescent/skyline.rs`](../../packages/treetime/src/coalescent/skyline.rs)                                                                                                  | complete      |
| [Relaxed Clock](timetree.md#relaxed-clock)                                      | well-known | timetree     | [`packages/treetime/src/timetree/optimization/relaxed_clock.rs`](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs)                                                                | complete      |
| [Convergence Monitoring](timetree.md#convergence-monitoring)                    | custom     | timetree     | [`packages/treetime/src/timetree/convergence/`](../../packages/treetime/src/timetree/convergence/)                                                                                                  | complete      |
| [Confidence Intervals](timetree.md#confidence-intervals)                        | well-known | timetree     | [`packages/treetime/src/timetree/confidence.rs`](../../packages/treetime/src/timetree/confidence.rs)                                                                                                | complete      |
| [FFT Convolution](distribution.md#fft-convolution)                              | well-known | distribution | [`packages/treetime-ops/src/convolution.rs`](../../packages/treetime-ops/src/convolution.rs)                                                                                                        | complete      |
| [Gaussian Convolution](distribution.md#gaussian-convolution)                    | well-known | distribution | [`packages/treetime-analytical/src/gaussian.rs`](../../packages/treetime-analytical/src/gaussian.rs)                                                                                                | complete      |
| [GTR Models](gtr.md#substitution-models)                                        | well-known | gtr          | [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)                                                                                                                | complete      |
| [Matrix Exponentiation](gtr.md#matrix-exponentiation)                           | well-known | gtr          | [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)                                                                                                                        | complete      |
| [Dense GTR Inference](gtr.md#dense-gtr-inference)                               | custom     | gtr          | [`packages/treetime/src/gtr/infer_gtr/common.rs`](../../packages/treetime/src/gtr/infer_gtr/common.rs)                                                                                              | complete      |
| [Parallel BFS](graph.md#parallel-bfs)                                           | well-known | graph        | [`packages/treetime-graph/src/breadth_first.rs`](../../packages/treetime-graph/src/breadth_first.rs)                                                                                                | complete      |
| [Newton-Raphson](optimization.md#newton-raphson-for-branch-length-optimization) | well-known | optimization | [`packages/treetime/src/optimize/dispatch.rs`](../../packages/treetime/src/optimize/dispatch.rs)                                                                                                    | complete      |
| [Poisson Indel Count](indel-models.md)                                          | custom     | indel        | [`packages/treetime/src/optimize/indel.rs`](../../packages/treetime/src/optimize/indel.rs)                                                                                                          | complete      |
| [TKF91](../reports/indel-models/3-single-residue.md)                            | well-known | indel        | -                                                                                                                                                                                                   | unimplemented |
| [PIP](../reports/indel-models/3-single-residue.md)                              | well-known | indel        | -                                                                                                                                                                                                   | unimplemented |

---

## File Index

| Domain       | Files                                                                                                                                                                                                                                                                                          | Algorithms                                                    |
| ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| Ancestral    | [`packages/treetime/src/ancestral/`](../../packages/treetime/src/ancestral/), [`packages/treetime/src/partition/marginal_*.rs`](../../packages/treetime/src/partition/)                                                                                                                        | Fitch, Marginal ML                                            |
| Clock        | [`packages/treetime/src/clock/`](../../packages/treetime/src/clock/)                                                                                                                                                                                                                           | WLS, regression, Brent, best root, outlier detection          |
| Timetree     | [`packages/treetime/src/commands/timetree/`](../../packages/treetime/src/commands/timetree/)                                                                                                                                                                                                   | Belief propagation, coalescent, convergence, relaxed clock    |
| Mugration    | [`packages/treetime/src/mugration/`](../../packages/treetime/src/mugration/), [`packages/treetime/src/partition/marginal_discrete.rs`](../../packages/treetime/src/partition/marginal_discrete.rs), [`packages/treetime/src/gtr/refinement.rs`](../../packages/treetime/src/gtr/refinement.rs) | Discrete marginal, GTR construction + refinement              |
| Distribution | [`packages/treetime-ops/src/`](../../packages/treetime-ops/src/), [`packages/treetime-analytical/src/`](../../packages/treetime-analytical/src/), [`packages/treetime-distribution/src/`](../../packages/treetime-distribution/src/)                                                           | Convolution, multiplication, analytical                       |
| GTR          | [`packages/treetime/src/gtr/`](../../packages/treetime/src/gtr/)                                                                                                                                                                                                                               | JC69, K80, HKY85, TN93, JTT92, inference (sparse + dense)     |
| Graph        | [`packages/treetime-graph/src/`](../../packages/treetime-graph/src/)                                                                                                                                                                                                                           | BFS, DFS, path finding                                        |
| Optimization | [`packages/treetime/src/optimize/`](../../packages/treetime/src/optimize/), [`packages/treetime-grid/src/`](../../packages/treetime-grid/src/)                                                                                                                                                 | Newton-Raphson, interpolation                                 |
| Indel Models | [`packages/treetime/src/optimize/indel.rs`](../../packages/treetime/src/optimize/indel.rs)                                                                                                                                                                                                     | Poisson indel count (implemented); TKF91, PIP (unimplemented) |

---

## References

### Primary

Sagulenko, P., Puller, V., & Neher, R.A. (2018). "TreeTime: Maximum-likelihood phylodynamic analysis." Virus Evolution, 4(1):vex042. https://doi.org/10.1093/ve/vex042

### Phylogenetics

- Felsenstein (1981). "Evolutionary trees from DNA sequences." J Mol Evol, 17(6):368-376. https://doi.org/10.1007/BF01734359
- Felsenstein (2003). "Inferring Phylogenies." Sinauer Associates. ISBN 978-0-87893-177-4
- Yang (2006). "Computational Molecular Evolution." Oxford University Press. ISBN 978-0-19-856702-8

### Substitution Models

- Jukes & Cantor (1969). "Evolution of Protein Molecules." In Munro (ed.), Mammalian Protein Metabolism, vol. 3, pp. 21-132. Academic Press
- Kimura (1980). "A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences." J Mol Evol, 16(2):111-120. https://doi.org/10.1007/BF01731581
- Hasegawa, Kishino & Yano (1985). "Dating of the human-ape splitting by a molecular clock of mitochondrial DNA." J Mol Evol, 22(2):160-174. https://doi.org/10.1007/BF02101694
- Tamura & Nei (1993). "Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees." Mol Biol Evol, 10(3):512-526. https://doi.org/10.1093/oxfordjournals.molbev.a040023
- Jones, Taylor & Thornton (1992). "The rapid generation of mutation data matrices from protein sequences." CABIOS, 8(3):275-282. https://doi.org/10.1093/bioinformatics/8.3.275

### Coalescent

- Kingman (1982). "The coalescent." Stochastic Processes and their Applications, 13(3):235-248. https://doi.org/10.1016/0304-4149(82)90011-4
- Strimmer & Pybus (2001). "Exploring the demographic history of DNA sequences using the generalized skyline plot." Mol Biol Evol, 18(12):2298-2305. https://doi.org/10.1093/oxfordjournals.molbev.a003776
- Drummond, Rambaut, Shapiro & Pybus (2005). "Bayesian coalescent inference of past population dynamics from molecular sequences." Mol Biol Evol, 22(5):1185-1192. https://doi.org/10.1093/molbev/msi103
- Minin, Bloomquist & Suchard (2008). "Smooth skyride through a rough skyline: Bayesian coalescent-based inference of population dynamics." Mol Biol Evol, 25(7):1459-1471. https://doi.org/10.1093/molbev/msn090

### Relaxed Clock

- Thorne, Kishino & Painter (1998). "Estimating the rate of evolution of the rate of molecular evolution." Mol Biol Evol, 15(12):1647-1657. https://doi.org/10.1093/oxfordjournals.molbev.a025892
- Drummond, Ho, Phillips & Rambaut (2006). "Relaxed phylogenetics and dating with confidence." PLOS Biology, 4(5):e88. https://doi.org/10.1371/journal.pbio.0040088
- Lepage, Bryant, Philippe & Lartillot (2007). "A general comparison of relaxed molecular clock models." Mol Biol Evol, 24(12):2669-2680. https://doi.org/10.1093/molbev/msm193

### Numerical Methods

- Brent (1973). "Algorithms for Minimization Without Derivatives." Prentice-Hall. ISBN 978-0-13-022335-7
- Pearl (1988). "Probabilistic Reasoning in Intelligent Systems: Networks of Plausible Inference." Morgan Kaufmann. ISBN 978-0-934613-73-2
- Nocedal & Wright (2006). "Numerical Optimization." 2nd ed. Springer. ISBN 978-0-387-30303-1
- Cormen, Leiserson, Rivest & Stein (2022). "Introduction to Algorithms." 4th ed. MIT Press. ISBN 978-0-262-04630-5

---

## Appendix: Classification

### Well-Known (Textbook/Published)

Fitch Parsimony, Felsenstein Pruning, Joint ML/Viterbi, Belief Propagation, Kingman Coalescent, Skyline Plot, Relaxed Clock DP, IQD Outlier Detection, Confidence Intervals, JC69/K80/F81/HKY85/T92/TN93, JTT92, Matrix Exponentiation, Newton-Raphson, Brent's Method, Golden Section Search, FFT Convolution, BFS/DFS, Piecewise Linear Interpolation

### Custom/TreeTime-Specific

Sufficient Statistics WLS, Sparse Marginal Reconstruction, Dense GTR Inference, Polytomy Resolution (Greedy), GTR Inference (EM-like), Lazy Normalization Multiplication, ScaledDistribution, Branch Point Cost Function, Variance Model, Best Root Search, Convergence Monitoring, Branch Length Distributions
