# Indel models in phylogenetics

> This document is AI-generated and may contain factual errors, misattributions, or outdated information. Verify claims against cited sources before relying on them.

A survey of insertion-deletion (indel) modeling approaches in phylogenetics: mathematical foundations, pair HMM and transducer formalisms, software implementations, and applicability to branch length optimization in TreeTime. Covers 12 models from gap-as-missing-data through full mechanistic birth-death processes. 30 verified references.

## Reading guide

Start with Chapter 1 for the problem statement and notation. Chapter 2 covers gap treatment approaches actually used in practice: missing data, the Poisson count model (implemented in TreeTime v1), and affine gap penalties. Chapters 3 and 4 cover the mechanistic birth-death models (TKF91, TKF92, RS07, PIP, GGI, SID) that provide theoretical foundations but are not implemented in TreeTime. Chapter 5 covers empirical rate models and parsimony approaches. Chapter 6 compares all models, analyzes their impact on branch lengths, and surveys 17 software implementations.

The [Glossary](glossary.md) defines every domain-specific term used across chapters. Complex terms are also defined inline on first use.

## Chapters

- [1. Introduction](1-introduction.md)
- [2. Gap treatment and counting models](2-gap-treatment.md)
- [3. Single-residue birth-death models](3-single-residue.md)
- [4. Multi-residue extensions and approximations](4-multi-residue.md)
- [5. Empirical models and parsimony](5-empirical.md)
- [6. Comparative analysis and software landscape](6-comparison.md)
- [Glossary](glossary.md)
- [References](references.md)
