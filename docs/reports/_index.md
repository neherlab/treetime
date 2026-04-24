# Reports

Technical reference books for TreeTime v1 algorithms and implementation. Each book is self-contained with its own table of contents, glossary, and references.

## Books

### [Iterative tree refinement](iterative-tree-refinement/_index.md)

How TreeTime refines phylogenetic trees: optimizing branch lengths, detecting and collapsing zero-length branches, resolving polytomies, and integrating these operations into convergent iteration loops. 10 chapters covering substitution models, tree likelihood, ancestral reconstruction, per-edge optimization, zero-length branches, polytomy resolution, initial estimation, the iteration loop, and implementation status. 74 verified references. Written for software engineers and bioinformatics students.

### [Command relationships: ancestral, prune, optimize, timetree](command-relationships/_index.md)

Scientific and architectural relationships between the four main tree-refinement commands. Covers current Venn diagram of shared components, the role of `ancestral` as the shared foundation, comparison of the two sibling EM-like loops (optimize vs timetree), current vs ideal architecture, and a gap table mapping six known issues to required changes.

### [Ancestral vs mugration comparison](ancestral-mugration-comparison/_index.md)

Side-by-side comparison of the `ancestral` and `mugration` commands in both v0 (Python) and v1 (Rust). Covers the shared Felsenstein pruning core, divergent pipelines (one-shot vs iterative GTR refinement), data representations (multi-position sequences vs single-position discrete traits), partition types, message-passing implementations, output formats, known issues, and architectural observations. Extends the command-relationships report by documenting mugration's position as a sibling outside the ancestral-prune-optimize-timetree hierarchy.

### [Indel models in phylogenetics](indel-models/_index.md)

Survey of insertion-deletion (indel) modeling approaches: mathematical foundations, pair HMM and transducer formalisms, software implementations, and applicability to branch length optimization. Covers 12 models from gap-as-missing-data through full mechanistic birth-death processes (TKF91, TKF92, PIP, GGI). Includes per-tool writeups of 17 phylogenetic software packages with implementation details. 6 chapters, glossary, 30 verified references.

### [Optimization methods](optimization-methods/_index.md)

Numerical optimization algorithms in TreeTime v1: per-edge branch length optimization (Newton-Raphson, Brent, cubic Hermite spline), coalescent skyline estimation (GMRF, NelderMead vs LBFGS), EM convergence, and the argmin Rust crate. Systematic comparison across RAxML-NG, IQ-TREE, PhyML, BEAST, and TreeTime v0/v1 with code locations. Includes a v1 optimization audit with 17 items inventoried, 7 inconsistencies identified, and 8 prioritized refactoring proposals. 7 chapters, glossary, 36 verified references.

### [Sparse substitution accessors](sparse-subs-accessors.md)

Production call sites for `fitch_subs` and `marginal_subs` on `SparseEdgePartition`. 22 fitch call sites across ancestral, marginal passes, sparse reroot, prune, topology cleanup, optimize, and GTR inference. 6 marginal call sites across marginal passes, sparse reroot, and output serialization.

### [Indel rate re-estimation in the optimize loop](optimize-indel-rate-reestimation.md)

Analysis of whether the global indel rate $\hat{\mu}$ should be re-estimated each iteration of the optimize loop or computed once and cached. Covers the ECM convergence guarantee, phylogenetic ML precedent (PhyML, RAxML-NG, IQ-TREE), the `gtr.mu` analogy, empirical 2-cycle evidence on sc2/2844, and the relationship to convergence proposal items P1-P4 and P9. 5 references.

### [Zero branch length optimization](zero-branch-length-optimization.md)

How TreeTime v1 decides when an edge has zero optimal branch length, handles boundary conditions, and collapses degenerate edges. Covers the derivative-sign shortcut for unimodal models (JC69, F81, binary), post-dispatch reconciliation for non-unimodal models (K80, HKY85, TN93, GTR), grid search fallback, input domain protection, and topology cleanup via edge collapse. Decision matrix for all 6 model/indel/validity combinations.
