# Chapter 1: Introduction

[Back to index](README.md) | Next: [Chapter 2: Branch length optimization](2-branch-length-optimization.md)

## Scope

This book covers the numerical optimization algorithms used throughout TreeTime v1 (Rust). It complements the [Iterative tree refinement](../iterative-tree-refinement/README.md) book, which describes the pipeline architecture. Here we focus on the algorithms themselves: what they optimize, how they converge, how they compare to implementations in other phylogenetic tools, and where the current v1 implementation can be improved.

## Classification of optimization problems

TreeTime solves several distinct optimization problems, each with different structure:

**1D scalar optimization (derivative-based)**

- Per-edge branch length: maximize log-likelihood `L(t) = sum_i log(sum_c k_{ic} exp(lambda_c t))` with analytical first and second derivatives. The innermost and most frequently called optimization in the pipeline.
- This is where Newton-Raphson, Brent, and cubic Hermite spline methods apply.

**1D scalar optimization (derivative-free)**

- Root-split position: minimize chi-squared along an edge, `x in [0, 1]`. No analytical derivatives available.
- Coalescent Tc: minimize negative log-likelihood in log-space, `log(Tc) in [-20, 2]`. Derivatives available but not implemented.
- Polytomy merge time: maximize likelihood gain, `t in [parent_time, child_min_time]`. Derivatives not implemented.
- These use Brent's method or golden section search via the argmin crate.

**Multi-dimensional optimization**

- Skyline coalescent: minimize penalized negative log-likelihood over n grid points of `log(Tc)`. Currently uses Nelder-Mead (derivative-free), but analytical gradients are tractable. 10-20 dimensions.

**Fixed-point iteration**

- GTR model parameters: alternating closed-form updates of rate matrix W, frequencies pi, and rate mu until convergence. Not a standard optimization problem.

**Analytical (closed-form)**

- Clock rate estimation: weighted least squares via sufficient statistics. Single-pass, no iteration.
- Relaxed clock gamma: quadratic penalty minimization via two-pass tree traversal. O(n), no iteration within a single step.

**Pipeline iteration (EM-like)**

- Outer optimize loop: alternates marginal reconstruction and branch length optimization with damping.
- Timetree refinement loop: iterates over time estimation, rerooting, branch optimization, ancestral reconstruction, and coalescent re-optimization.

## Relationship to the iterative tree refinement book

The [iterative tree refinement book](../iterative-tree-refinement/README.md) describes what each optimization does in the pipeline and how they interact. This book describes the algorithms themselves: convergence properties, computational cost, failure modes, and how competing tools solve the same problems. The two books share a reference list with significant overlap.

Cross-references between the books use the pattern `../iterative-tree-refinement/<chapter>.md`.

## The argmin crate

TreeTime v1 uses the [argmin](https://docs.rs/argmin/0.10.0/argmin/) Rust crate as its optimization framework. argmin provides an extensible `Solver` trait, an `Executor` for running solvers with observers and convergence criteria, and built-in solvers including `BrentOpt`, `BrentRoot`, `GoldenSectionSearch`, `NelderMead`, `LBFGS`, and `Newton`. Chapter 6 covers argmin's architecture and how to write custom solvers.

## Convention

Throughout this book, optimization items from the [audit](7-audit.md) are referenced by their identifiers: E1-E6 for existing argmin usage, H1-H11 for hand-rolled implementations, I1-I7 for inconsistencies, and P1-P8 for refactoring proposals.
