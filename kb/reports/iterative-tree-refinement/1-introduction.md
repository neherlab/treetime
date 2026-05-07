# Chapter 1: Introduction

[Back to index](README.md) | Next: [Chapter 2: Substitution models](2-substitution-models.md)

## What is a phylogenetic tree?

A phylogenetic tree is a diagram showing how organisms are related through evolution. In molecular phylogenetics, organisms are represented by their DNA (or protein) sequences. The tree has three types of elements:

- **Leaves** (tips): observed sequences. Each leaf is a sample -- a virus isolate, a bacterial strain, a patient's tumor, or any other sequenced entity.
- **Internal nodes**: hypothetical ancestors whose sequences are not directly observed but can be inferred.
- **Branches** (edges): connections between nodes. Each branch represents an evolutionary lineage -- a line of descent from ancestor to descendant.

```
        root
       /    \
      /      \
    AB        CD
   / \       / \
  A   B     C   D
```

A, B, C, D are observed sequences (leaves). AB and CD are inferred ancestors (internal nodes). The root is the common ancestor of all samples.

## What are branch lengths?

Each branch has a **length** measuring the amount of evolutionary change. Branch lengths are measured in **expected substitutions per site** -- the average number of nucleotide changes per alignment position.

A branch length of 0.01 means roughly 1 substitution per 100 positions. A branch length of 0 means the parent and child are identical -- no observed evolution.

Branch lengths encode three things:

- **How much evolution happened**: longer branches mean more genetic change
- **How much time passed** (combined with a molecular clock): branch length / mutation rate = time
- **How confident we are in the branching order**: very short internal branches mean the topology is uncertain at that split

## What is tree refinement?

Tree-building programs (IQ-TREE, RAxML, FastTree) produce trees with approximate branch lengths and fixed binary topology. Tree refinement improves the tree by adjusting two things:

**Branch lengths** (edge weights): the optimizer finds the lengths that make the observed sequences most probable under the substitution model. This is a continuous optimization problem -- adjusting numbers on a fixed structure.

**Tree topology** (edge structure): the optimizer detects and removes unsupported internal branches (creating polytomies) and resolves polytomies where the data supports binary structure (adding internal nodes). This is a discrete topology-modification problem -- adding and removing nodes.

These two operations are interdependent:

1. Optimizing branch lengths reveals that some branches should be zero
2. Collapsing zero-length branches changes the topology (creates polytomies)
3. Resolving polytomies introduces new edges with unknown lengths
4. Those new edges need branch length optimization
5. Repeat until both stabilize

This alternation is the **iteration loop** described in [Chapter 9](9-iteration-loop.md).

## The topology problem

Tree builders produce **fully bifurcating** trees: every internal node has exactly two children. This is a modeling assumption, not a biological fact. When four virus sequences diverged simultaneously during an outbreak, the true history is a four-way split (a **polytomy**), but the tree builder picks one of 15 possible binary resolutions:

```
True history (polytomy):          Tree builder output (arbitrary):

      root                              root
    / | | \                            /    \
   A  B  C  D                       AB      CD
                                    / \    / \
                                   A   B  C   D
```

The internal branches root-AB and root-CD have near-zero length because they represent structure the data does not support. Two operations address this:

1. **Prune zero-length branches**: collapse unsupported internal edges, creating polytomies. This admits that the data cannot resolve a particular split. See [Chapter 6](6-zero-length-branches.md).

2. **Resolve polytomies by shared mutations**: scan children of a polytomy for shared substitutions and introduce new internal nodes grouping children that share mutations. This uses sequence signal to create binary splits where the data supports them. See [Chapter 7](7-polytomy-resolution.md).

These are complementary: pruning removes unsupported structure, and shared-mutation merging adds supported structure.

## Where TreeTime fits

TreeTime is a phylodynamic analysis tool used in the Nextstrain platform for real-time pathogen surveillance. Its primary users are the Nextstrain pathogen pipelines (flu, SARS-CoV-2, Ebola, Zika, and others) via the Augur toolkit. TreeTime was introduced by Sagulenko, Puller, and Neher (2018).

TreeTime has three commands relevant to tree refinement:

- **`optimize`**: branch length optimization. Takes a tree and alignment, refines branch lengths by maximum likelihood. Topology cleanup (pruning and polytomy resolution) should be integrated here but is not yet fully implemented.
- **`prune`**: standalone tree cleanup. Removes short branches, empty branches, or named nodes. Can merge sibling branches sharing mutations.
- **`timetree`**: the main phylodynamic pipeline. Combines branch length optimization with molecular clock inference, polytomy resolution, coalescent modeling, and ancestral reconstruction in an iterative loop.

## What this book covers

The remaining chapters explain each component of tree refinement:

- [Chapter 2](2-substitution-models.md): How DNA substitution is modeled mathematically
- [Chapter 3](3-tree-likelihood.md): How the likelihood of a tree is computed
- [Chapter 4](4-ancestral-reconstruction.md): How ancestral sequences are inferred
- [Chapter 5](5-branch-length-optimization.md): How individual branch lengths are optimized
- [Chapter 6](6-zero-length-branches.md): When and how zero-length branches are detected and removed
- [Chapter 7](7-polytomy-resolution.md): How polytomies are resolved
- [Chapter 8](8-initial-estimation.md): How branch lengths are initialized before optimization
- [Chapter 9](9-iteration-loop.md): How these components fit together in the iteration loop
- [Chapter 10](10-status-and-recommendations.md): Current implementation status and recommendations

Each chapter begins with an accessible introduction, then proceeds to mathematical detail, v0/v1 code references, and implementation advice.

## References

- Sagulenko, P., V. Puller, and R. A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042
