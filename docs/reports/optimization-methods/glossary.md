# Glossary

[Back to index](_index.md)

Terms used throughout this book. For phylogenetic terms not listed here, see the [Iterative tree refinement glossary](../iterative-tree-refinement/glossary.md).

## argmin

A Rust crate for numerical optimization providing the `Solver` trait, `Executor` runner, and built-in solvers (BrentOpt, BrentRoot, GoldenSectionSearch, NelderMead, LBFGS, Newton). TreeTime v1's optimization framework. https://docs.rs/argmin/0.10.0/argmin/

## Brent's method

A derivative-free 1D optimization algorithm combining golden-section search with successive parabolic interpolation. Guaranteed to converge within a bracket. Convergence order ~1.325 (tribonacci constant). Used by TreeTime v0 for branch lengths and by v1 (via argmin BrentOpt) for root-split, Tc, and polytomy merge-time optimization. Brent 1973.

## BrentOpt / BrentRoot

argmin's implementations of Brent's method. BrentOpt minimizes a function on a bounded interval. BrentRoot finds a zero crossing on a bounded interval. Both need only `CostFunction` (no derivatives).

## Convergence order

The rate at which an iterative method approaches the solution. Order 1 (linear): constant ratio of errors. Order ~1.325 (superlinear): Brent's parabolic interpolation. Order ~1.618 (superlinear): secant method. Order 2 (quadratic): Newton-Raphson. Order 3 (cubic): PhyML's Hermite spline. Higher order means fewer iterations but often higher per-iteration cost.

## Cost function

In argmin, a type implementing the `CostFunction` trait with `cost(&self, param: &P) -> Result<F, Error>`. Maps a parameter value to a scalar cost. The optimizer minimizes this.

## Cubic Hermite spline

PhyML's branch length optimization method. Given function values and first derivatives at two bracketing points, constructs a cubic polynomial and solves for its stationary points. Achieves cubic convergence using only first derivatives (no second derivatives needed).

## ECM (Expectation-Conditional Maximization)

An extension of EM where the M-step is replaced by sequential conditional maximizations over parameter subsets. Preserves the monotone likelihood guarantee. TreeTime's alternating reconstruction + branch optimization + damping follows this pattern. Meng and Rubin 1993.

## EM (Expectation-Maximization)

An iterative algorithm for maximum likelihood estimation with latent variables. E-step: compute posterior of latent variables. M-step: maximize expected complete-data log-likelihood. Guarantees monotone likelihood increase. Linear convergence rate governed by fraction of missing information. Dempster, Laird, and Rubin 1977.

## Executor

argmin's runner that wraps a `Solver` and a problem, manages iteration, convergence criteria, observers, and checkpointing. Configured via `Executor::new(problem, solver).configure(|state| state.max_iters(100)).run()`.

## GMRF (Gaussian Markov Random Field)

A multivariate Gaussian distribution with a sparse (tridiagonal) precision matrix Q. Used as a smoothness prior on log-transformed population size trajectories in coalescent skyline models. The penalty `theta^T Q theta` penalizes differences between adjacent grid values. Minin et al. 2008.

## Golden section search

A derivative-free bracketed optimization method that narrows the bracket by the golden ratio (0.618) each iteration. Linear convergence. Guaranteed to find the minimum of a unimodal function within the bracket. Used by TreeTime v1 (via argmin GoldenSectionSearch) for root-split optimization.

## LBFGS (Limited-memory BFGS)

A quasi-Newton optimization method that approximates the Hessian inverse using m recent gradient pairs. Memory O(mn) instead of O(n^2) for full BFGS. Requires gradient but not Hessian. Superlinear convergence for smooth functions. Liu and Nocedal 1989.

## Nelder-Mead

A derivative-free optimization method using a simplex of n+1 vertices in n dimensions. Operations: reflection, expansion, contraction, shrink. No convergence guarantee above 2D (McKinnon 1998). Currently used by TreeTime v1 for skyline optimization. Nelder and Mead 1965.

## Newton-Raphson

A root-finding/optimization method using first and second derivatives: `x_new = x - f'(x)/f''(x)`. Quadratic convergence near the optimum. Requires the Hessian to be negative definite (function locally concave). The standard method for phylogenetic branch length optimization (RAxML-NG, IQ-TREE, TreeTime v1).

## Observer

In argmin, a type implementing the `Observe` trait that receives notifications during optimization. Used for logging, progress tracking, and diagnostics. Attached via `Executor.add_observer(obs, ObserverMode::Every(10))`.

## SLSQP (Sequential Least Squares Programming)

A constrained optimization algorithm that solves a sequence of QP subproblems using least-squares. Handles equality, inequality, and bound constraints. Superlinear convergence. Used by TreeTime v0 (via scipy) for skyline optimization. Kraft 1988.

## Solver trait

argmin's core abstraction. Implementing `Solver<O, I>` requires `next_iter(&mut self, problem, state) -> (state, kv)` and optionally `init`, `terminate`. Enables custom optimization algorithms to use argmin's executor infrastructure.

## Unimodality

A function with at most one local maximum. The per-edge phylogenetic log-likelihood is provably unimodal under JC69 and F81 models, but can have multiple maxima under K2P and more complex models. Dinh and Matsen 2017.
