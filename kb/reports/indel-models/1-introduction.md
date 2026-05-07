# Chapter 1: Introduction

[Back to index](README.md) | Next: [Chapter 2: Gap treatment and counting models](2-gap-treatment.md)

## The problem

Standard phylogenetic substitution models (JC69 through GTR) treat sequence length as fixed and model only character-state changes at individual sites. Insertions and deletions (indels) are handled separately: either discarded as missing data, penalized heuristically during alignment, or modeled by dedicated stochastic processes.

This creates a blind spot in branch length optimization: a branch with zero substitutions but one or more indels is assigned zero length, collapsing topology that the indel evidence supports. TreeTime v1 addresses this with a Poisson indel count model (Chapter 2), but the broader field offers a spectrum of approaches with different tradeoffs.

## Scope

This report surveys indel modeling approaches in phylogenetics, organized from simplest to most complex:

- Gap penalties and counting models (Chapter 2): gaps as missing data, Poisson count, affine gap penalty
- Single-residue birth-death processes (Chapter 3): TKF91, PIP
- Multi-residue extensions (Chapter 4): TKF92, RS07, GGI, SID
- Empirical rate models and parsimony (Chapter 5): SIM/RIM, indelMaP, binary encoding
- Comparative analysis (Chapter 6): model hierarchy, impact on branch lengths, software landscape

Each model section covers mathematical formulation, computational complexity, software implementations, and applicability to TreeTime's branch length optimization.

## Notation

| Symbol     | Meaning                                                       |
| :--------- | :------------------------------------------------------------ |
| $\lambda$  | Insertion (birth) rate per site per unit time                 |
| $\mu$      | Deletion (death) rate per site per unit time                  |
| $t$        | Branch length (evolutionary distance, substitutions per site) |
| $k$        | Number of observed indel events on an edge                    |
| $L$        | Alignment length                                              |
| $N$        | Number of taxa (sequences)                                    |
| $\beta(t)$ | TKF91 absorption probability (see Chapter 3)                  |
