# Reports

Technical reference books for TreeTime v1 algorithms and implementation. Each book is self-contained with its own table of contents, glossary, and references.

## Books

### [Iterative tree refinement](iterative-tree-refinement/_index.md)

How TreeTime refines phylogenetic trees: optimizing branch lengths, detecting and collapsing zero-length branches, resolving polytomies, and integrating these operations into convergent iteration loops. 10 chapters covering substitution models, tree likelihood, ancestral reconstruction, per-edge optimization, zero-length branches, polytomy resolution, initial estimation, the iteration loop, and implementation status. 74 verified references. Written for software engineers and bioinformatics students.

### [Command relationships: prune, optimize, timetree](command-relationships/_index.md)

Scientific and architectural relationships between the three main tree-refinement commands. Covers current Venn diagram of shared components, comparison of the two sibling EM-like loops (optimize vs timetree), current vs ideal architecture, and a gap table mapping six known issues to required changes.

### [Optimization methods](optimization-methods/_index.md)

Numerical optimization algorithms in TreeTime v1: per-edge branch length optimization (Newton-Raphson, Brent, cubic Hermite spline), coalescent skyline estimation (GMRF, NelderMead vs LBFGS), EM convergence, and the argmin Rust crate. Systematic comparison across RAxML-NG, IQ-TREE, PhyML, BEAST, and TreeTime v0/v1 with code locations. Includes a v1 optimization audit with 17 items inventoried, 7 inconsistencies identified, and 8 prioritized refactoring proposals. 7 chapters, glossary, 36 verified references.
