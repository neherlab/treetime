# Algorithm Inventory - TreeTime v1 (Rust) and v0 (Python)

> **Valid for commit**: `41c0fb16` (rust branch, 2026-03-10)

> **Note**: This document was AI-generated via parallel agent analysis of the codebase. Verify critical details against source code.

Accepted-scope rule: this inventory documents algorithms and algorithmic workflows that are implemented in v1 or explicitly tracked as unimplemented port targets. Candidate extensions, disputed classifications, and not-yet-accepted behavioral changes belong in [`docs/port-proposals/`](../port-proposals/_index.md), not here.

## Summary

| Domain                                      | v1      | Unimplemented | Well-Known | Custom  |
| ------------------------------------------- | ------- | ------------- | ---------- | ------- |
| [Ancestral Reconstruction](ancestral.md)    | 2       | 1             | 3          | 0       |
| [Clock Inference](clock.md)                 | 14      | 4             | 8          | 8       |
| [Timetree Inference](timetree.md)           | 20      | 1             | 11         | 6       |
| [Mugration](mugration.md)                   | 3       | 1             | 2          | 0       |
| [Distribution/Convolution](distribution.md) | 17      | 4             | 12         | 5       |
| [GTR Substitution Models](gtr.md)           | 12      | 3             | 8          | 3       |
| [Graph Traversal](graph.md)                 | 4       | 0             | 4          | 0       |
| [Numerical Optimization](optimization.md)   | 4       | 0             | 5          | 0       |
| **Total**                                   | **~76** | **14**        | **~53**    | **~22** |

Full unimplemented algorithm details: [unimplemented.md](unimplemented.md)

---

## Quick Reference

| Algorithm                                                          | Type       | Domain       | Location                                                                                                                                               | Status        |
| ------------------------------------------------------------------ | ---------- | ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------- |
| [Fitch Parsimony](ancestral.md#fitch-parsimony)                    | well-known | ancestral    | [`packages/treetime/src/commands/ancestral/fitch.rs`](../../packages/treetime/src/commands/ancestral/fitch.rs)                                         | complete      |
| [Marginal ML](ancestral.md#marginal-ml)                            | well-known | ancestral    | [`packages/treetime/src/representation/partition/marginal_*.rs`](../../packages/treetime/src/representation/partition/)                                | complete      |
| [Discrete Marginal](mugration.md#discrete-marginal-reconstruction) | well-known | mugration    | [`packages/treetime/src/commands/mugration/discrete_marginal.rs`](../../packages/treetime/src/commands/mugration/discrete_marginal.rs)                 | complete      |
| [Joint ML](unimplemented.md#joint-ml)                              | well-known | ancestral    | -                                                                                                                                                      | unimplemented |
| [WLS Sufficient Stats](clock.md#wls-sufficient-stats)              | custom     | clock        | [`packages/treetime/src/commands/clock/clock_set.rs`](../../packages/treetime/src/commands/clock/clock_set.rs)                                         | complete      |
| [Tree Regression](clock.md#tree-regression)                        | well-known | clock        | [`packages/treetime/src/commands/clock/clock_regression.rs`](../../packages/treetime/src/commands/clock/clock_regression.rs)                           | complete      |
| [Brent's Method](clock.md#brents-method)                           | well-known | clock        | [`packages/treetime/src/commands/clock/find_best_root/method_brent.rs`](../../packages/treetime/src/commands/clock/find_best_root/method_brent.rs)     | complete      |
| [IQD Outlier Detection](clock.md#iqd-outlier-detection)            | well-known | clock        | [`packages/treetime/src/commands/clock/clock_filter.rs`](../../packages/treetime/src/commands/clock/clock_filter.rs)                                   | complete      |
| [Belief Propagation](timetree.md#belief-propagation)               | well-known | timetree     | [`packages/treetime/src/commands/timetree/inference/`](../../packages/treetime/src/commands/timetree/inference/)                                       | complete      |
| [Kingman Coalescent](timetree.md#kingman-coalescent)               | well-known | timetree     | [`packages/treetime/src/commands/timetree/coalescent/`](../../packages/treetime/src/commands/timetree/coalescent/)                                     | complete      |
| [Skyline Coalescent](timetree.md#skyline-coalescent)               | well-known | timetree     | [`packages/treetime/src/commands/timetree/coalescent/skyline.rs`](../../packages/treetime/src/commands/timetree/coalescent/skyline.rs)                 | complete      |
| [Relaxed Clock](timetree.md#relaxed-clock)                         | well-known | timetree     | [`packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs`](../../packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs) | complete      |
| [Convergence Monitoring](timetree.md#convergence-monitoring)       | custom     | timetree     | [`packages/treetime/src/commands/timetree/convergence/`](../../packages/treetime/src/commands/timetree/convergence/)                                   | complete      |
| [Confidence Intervals](timetree.md#confidence-intervals)           | well-known | timetree     | [`packages/treetime/src/commands/timetree/output/confidence.rs`](../../packages/treetime/src/commands/timetree/output/confidence.rs)                   | complete      |
| [FFT Convolution](distribution.md#fft-convolution)                 | well-known | distribution | [`packages/treetime-ops/src/convolution.rs`](../../packages/treetime-ops/src/convolution.rs)                                                           | complete      |
| [Gaussian Convolution](distribution.md#gaussian-convolution)       | well-known | distribution | [`packages/treetime-analytical/src/gaussian.rs`](../../packages/treetime-analytical/src/gaussian.rs)                                                   | complete      |
| [GTR Models](gtr.md#substitution-models)                           | well-known | gtr          | [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)                                                                   | complete      |
| [Matrix Exponentiation](gtr.md#matrix-exponentiation)              | well-known | gtr          | [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)                                                                           | complete      |
| [Dense GTR Inference](gtr.md#dense-gtr-inference)                  | custom     | gtr          | [`packages/treetime/src/gtr/infer_gtr/dense.rs`](../../packages/treetime/src/gtr/infer_gtr/dense.rs)                                                   | complete      |
| [Parallel BFS](graph.md#parallel-bfs)                              | well-known | graph        | [`packages/treetime-graph/src/breadth_first.rs`](../../packages/treetime-graph/src/breadth_first.rs)                                                   | complete      |
| [Newton-Raphson](optimization.md#newton-raphson)                   | well-known | optimization | [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)                     | complete      |

---

## File Index

| Domain       | Files                                                                                                                                                                                                                                | Algorithms                                                 |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------- |
| Ancestral    | [`packages/treetime/src/commands/ancestral/`](../../packages/treetime/src/commands/ancestral/), [`packages/treetime/src/representation/partition/marginal_*.rs`](../../packages/treetime/src/representation/partition/)              | Fitch, Marginal ML                                         |
| Clock        | [`packages/treetime/src/commands/clock/`](../../packages/treetime/src/commands/clock/)                                                                                                                                               | WLS, regression, Brent, best root, outlier detection       |
| Timetree     | [`packages/treetime/src/commands/timetree/`](../../packages/treetime/src/commands/timetree/)                                                                                                                                         | Belief propagation, coalescent, convergence, relaxed clock |
| Mugration    | [`packages/treetime/src/commands/mugration/`](../../packages/treetime/src/commands/mugration/), [`packages/treetime/src/representation/partition/discrete.rs`](../../packages/treetime/src/representation/partition/discrete.rs)     | Discrete marginal, GTR construction                        |
| Distribution | [`packages/treetime-ops/src/`](../../packages/treetime-ops/src/), [`packages/treetime-analytical/src/`](../../packages/treetime-analytical/src/), [`packages/treetime-distribution/src/`](../../packages/treetime-distribution/src/) | Convolution, multiplication, analytical                    |
| GTR          | [`packages/treetime/src/gtr/`](../../packages/treetime/src/gtr/)                                                                                                                                                                     | JC69, K80, HKY85, TN93, JTT92, inference (sparse + dense)  |
| Graph        | [`packages/treetime-graph/src/`](../../packages/treetime-graph/src/)                                                                                                                                                                 | BFS, DFS, path finding                                     |
| Optimization | [`packages/treetime/src/commands/optimize/`](../../packages/treetime/src/commands/optimize/), [`packages/treetime-grid/src/`](../../packages/treetime-grid/src/)                                                                     | Newton-Raphson, interpolation                              |

---

## References

### Primary

Sagulenko, P., Puller, V., & Neher, R.A. (2018). "TreeTime: Maximum-likelihood phylodynamic analysis." Virus Evolution, 4(1):vex042. https://doi.org/10.1093/ve/vex042

### Phylogenetics

- Felsenstein (1981). "Evolutionary trees from DNA sequences." J Mol Evol, 17(6):368-376
- Felsenstein (2003). "Inferring Phylogenies." Sinauer Associates
- Yang (2006). "Computational Molecular Evolution." Oxford University Press

### Substitution Models

- Jukes & Cantor (1969). Evolution of Protein Molecules
- Kimura (1980). J Mol Evol 16(2):111-120
- Hasegawa et al. (1985). J Mol Evol 22(2):160-174
- Tamura & Nei (1993). Mol Biol Evol 10(3):512-526
- Jones et al. (1992). CABIOS 8(3):275-282

### Coalescent

- Kingman (1982). "The coalescent." Stochastic Processes and Applications, 13(3):235-248
- Strimmer & Pybus (2001). Mol Biol Evol 18(12):2298-2305
- Drummond et al. (2005). "Bayesian coalescent inference." Mol Biol Evol 22(5):1185-1192
- Minin et al. (2008). Mol Biol Evol 25(7):1459-1471

### Relaxed Clock

- Thorne, Kishino & Painter (1998). "Estimating the rate of evolution of the rate of molecular evolution." Mol Biol Evol 15(12):1647-1657
- Drummond et al. (2006). "Relaxed phylogenetics and dating with confidence." PLOS Biology 4(5):e88
- Lepage et al. (2007). "A general comparison of relaxed molecular clock models." Mol Biol Evol 24(12):2669-2680

### Numerical Methods

- Brent (1973). "Algorithms for Minimization Without Derivatives." Prentice-Hall
- Pearl (1988). "Probabilistic Reasoning in Intelligent Systems." Morgan Kaufmann
- Nocedal & Wright. "Numerical Optimization." Springer
- CLRS. "Introduction to Algorithms." MIT Press

---

## Appendix: Classification

### Well-Known (Textbook/Published)

Fitch Parsimony, Felsenstein Pruning, Joint ML/Viterbi, Belief Propagation, Kingman Coalescent, Skyline Plot, Relaxed Clock DP, IQD Outlier Detection, Confidence Intervals, JC69/K80/F81/HKY85/T92/TN93, JTT92, Matrix Exponentiation, Newton-Raphson, Brent's Method, Golden Section Search, FFT Convolution, BFS/DFS, Piecewise Linear Interpolation

### Custom/TreeTime-Specific

Sufficient Statistics WLS, Sparse Marginal Reconstruction, Dense GTR Inference, Polytomy Resolution (Greedy), GTR Inference (EM-like), Lazy Normalization Multiplication, ScaledDistribution, Branch Point Cost Function, Variance Model, Best Root Search, Convergence Monitoring, Branch Length Distributions
