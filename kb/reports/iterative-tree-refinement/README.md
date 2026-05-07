# Iterative tree refinement in TreeTime

A reference on how TreeTime refines phylogenetic trees: optimizing branch lengths, detecting and collapsing zero-length branches, resolving polytomies, and integrating these operations into convergent iteration loops. Written for software engineers and bioinformatics students. All claims are supported by references to the scientific literature.

## Reading guide

Start with Chapter 1 for the biological motivation and the big picture. Chapters 2-4 build the mathematical foundation: substitution models, tree likelihood, and ancestral reconstruction. Chapters 5-8 cover the four operations that refine a tree: per-edge branch length optimization, zero-length branch detection and pruning, polytomy resolution, and initial branch length estimation. Chapter 9 explains how these operations fit together in the iteration loop, including damping, convergence theory, and the EM algorithm. Chapter 10 compares the current v0 and v1 implementations and lists open work.

The [Glossary](glossary.md) defines every technical term. Complex terms are also defined inline on first use in each chapter.

## Chapters

- [1. Introduction: trees, branch lengths, and why refinement matters](1-introduction.md)
- [2. Substitution models and sequence evolution](2-substitution-models.md)
- [3. Tree likelihood and the pruning algorithm](3-tree-likelihood.md)
- [4. Ancestral sequence reconstruction](4-ancestral-reconstruction.md)
- [5. Per-edge branch length optimization](5-branch-length-optimization.md)
- [6. Zero-length branches: detection and pruning](6-zero-length-branches.md)
- [7. Polytomy resolution](7-polytomy-resolution.md)
- [8. Initial branch length estimation](8-initial-estimation.md)
- [9. The iteration loop: putting it all together](9-iteration-loop.md)
- [10. Implementation status and recommendations](10-status-and-recommendations.md)
- [Glossary](glossary.md)
- [References](references.md)
