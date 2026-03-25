# Optimization methods in TreeTime

A reference on numerical optimization methods used throughout TreeTime v1: per-edge branch length optimization, coalescent skyline estimation, root-split search, polytomy resolution, and their integration via the argmin crate. Covers foundational algorithms, state-of-the-art approaches, implementations in competing phylogenetic tools (RAxML-NG, IQ-TREE, PhyML, BEAST), and actionable refactoring proposals for TreeTime v1. All claims are supported by references to the scientific literature or source code.

Companion to [Iterative tree refinement](../iterative-tree-refinement/_index.md), which covers the pipeline architecture. This book focuses on the optimization algorithms themselves.

## Reading guide

Start with Chapter 1 for the scope and the classification of optimization problems in TreeTime. Chapter 2 is the core: per-edge branch length optimization, comparing Newton-Raphson, Brent, and cubic Hermite spline across five tools. Chapter 3 covers coalescent skyline estimation with GMRF priors. Chapter 4 explains the outer EM-like loop and convergence theory. Chapter 5 covers the remaining optimization sub-problems (root-split, Tc, polytomy, HPD, GTR, relaxed clock). Chapter 6 documents the argmin Rust crate and how to write custom solvers. Chapter 7 is the audit: current inventory, inconsistencies across commands, and prioritized refactoring proposals.

## Chapters

- [1. Introduction: optimization in phylogenetics](1-introduction.md)
- [2. Branch length optimization](2-branch-length-optimization.md)
- [3. Coalescent skyline optimization](3-coalescent-skyline.md)
- [4. Outer loops and EM convergence](4-outer-loops.md)
- [5. Supporting optimizations](5-supporting-optimizations.md)
- [6. The argmin crate](6-argmin-crate.md)
- [7. Audit: inventory, inconsistencies, and proposals](7-audit.md)
- [Glossary](glossary.md)
- [References](references.md)
