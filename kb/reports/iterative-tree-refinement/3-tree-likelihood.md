# Chapter 3: Tree likelihood and the pruning algorithm

[Back to index](README.md) | Previous: [Chapter 2: Substitution models](2-substitution-models.md) | Next: [Chapter 4: Ancestral reconstruction](4-ancestral-reconstruction.md)

## What is tree likelihood?

The **likelihood** of a phylogenetic tree is the probability of observing the alignment data given the tree topology, branch lengths, and substitution model. Formally:

```
L(T, t, Q | D) = P(D | T, t, Q)
```

where T is the topology, t is the vector of branch lengths, Q is the substitution model, and D is the alignment. The tree with the highest likelihood is the maximum likelihood (ML) tree -- the tree that best explains the data. The ML approach to phylogenetics was introduced by Felsenstein (1973) and formalized with the pruning algorithm in Felsenstein (1981).

Computing this probability requires summing over all possible ancestral sequences at internal nodes, since those sequences are not observed. For a tree with n internal nodes, each with L alignment positions and s possible states per position, the naive computation requires `s^(n*L)` terms -- astronomically many.

Felsenstein's pruning algorithm (Felsenstein 1981) makes this tractable by exploiting the tree structure.

## The pruning algorithm

The pruning algorithm (also called the peeling algorithm) is a dynamic programming method that computes the likelihood in O(n \* s^2) time per site, where n is the number of nodes and s is the number of states (4 for nucleotides).

The key insight: the tree factorizes the probability. Each internal node's contribution depends only on its children, not on the rest of the tree. This means the sum over ancestral states can be computed node by node, from leaves to root.

### Setup

At each node u and each alignment position i, we maintain a **conditional likelihood vector** `L_u(i)` of length s. Entry `L_u(i)[a]` is the probability of observing all the data in u's subtree, given that node u has state a at position i.

### Leaf initialization

For a leaf node u with observed state `x_i` at position i:

```
L_u(i)[a] = 1  if a == x_i
L_u(i)[a] = 0  otherwise
```

If the leaf has an ambiguous state (e.g., N = any nucleotide), all entries are 1.

### Internal node recursion

For an internal node u with children v and w, connected by branches of length `t_v` and `t_w`:

```
L_u(i)[a] = [sum_b P(t_v)_{ab} * L_v(i)[b]] * [sum_c P(t_w)_{ac} * L_w(i)[c]]
```

Each child contributes a factor: the sum over all possible child states of (transition probability from parent state a to child state) times (child's conditional likelihood for that state). The two children's contributions are multiplied because they are conditionally independent given the parent state.

### Root likelihood

At the root, the site likelihood is:

```
L_site(i) = sum_a pi_a * L_root(i)[a]
```

where `pi_a` is the equilibrium frequency of state a. The total log-likelihood is:

```
log L = sum_i log L_site(i)
```

### Complexity

- Per site: O(n \* s^2) -- one matrix-vector product per node, n nodes
- Total: `O(n * s^2 * L)` -- L alignment positions
- For nucleotides (s=4): the s^2 = 16 factor is small, so the algorithm is effectively O(n \* L)

### Messages and profiles

In TreeTime's terminology, the conditional likelihood vectors are called **messages**. The backward pass (leaves to root) computes `msg_to_parent` at each node. The forward pass (root to leaves, described in [Chapter 4](4-ancestral-reconstruction.md)) computes `msg_to_child`. Together, these messages enable marginal reconstruction and per-edge branch length optimization.

The product of the two messages at a node gives the **posterior profile** -- the probability distribution over states at that node given all the data in the tree:

```
profile(u, i)[a] = msg_to_parent(u, i)[a] * msg_to_child(u, i)[a] / normalization
```

v1 code:

- Backward pass: `process_node_backward()` in [`packages/treetime/src/partition/marginal_passes.rs`](../../../packages/treetime/src/partition/marginal_passes.rs)
- Forward pass: `process_node_forward()` in the same file
- The `update_marginal()` function in [`packages/treetime/src/ancestral/marginal.rs`](../../../packages/treetime/src/ancestral/marginal.rs) orchestrates both passes

v0 code:

- `_ml_anc_joint()` and `_ml_anc_marginal()` in [`packages/legacy/treetime/treetime/treeanc.py`](../../../packages/legacy/treetime/treetime/treeanc.py)

## Numerical considerations

### Log-space computation

For long alignments, the product of per-site likelihoods underflows to zero in floating-point arithmetic. The standard solution is to work in log space:

```
log L = sum_i log L_site(i)
```

Each site likelihood is computed normally (it is a sum of s terms, not a product), then its log is taken. The total log-likelihood is a sum of logs.

Within a single site, the conditional likelihood vector entries can also underflow for deep trees with many nodes. The remedy is **scaling**: multiply each node's conditional likelihood vector by a scaling factor to keep values in a representable range, then account for the scaling factors in the final log-likelihood. This technique was introduced by Felsenstein (1981) and is standard in all ML phylogenetic software.

v1 code: The sparse representation in v1 uses log-space arithmetic with `logsumexp` normalization in [`packages/treetime/src/partition/marginal_helpers.rs`](../../../packages/treetime/src/partition/marginal_helpers.rs). The `softmax_with_log_norm()` function handles numerical stability for the backward and forward passes.

### Sparse vs dense representation

TreeTime supports two representations of sequence data on the tree:

**Dense**: stores the full probability vector at every position for every node. `O(n * L * s)` memory. Required when branches are long (many positions differ from the reference) or when the posterior profiles are needed for downstream analysis.

**Sparse**: stores only positions that differ from a reference sequence. Invariant positions (same state across all taxa) carry no information for branch length optimization and can be skipped. `O(n * V * s)` memory, where V is the number of variable positions (V << L for closely related sequences like viral isolates).

v1 code:

- Dense: `PartitionMarginalDense` in [`packages/treetime/src/partition/marginal_dense.rs`](../../../packages/treetime/src/partition/marginal_dense.rs)
- Sparse: `PartitionMarginalSparse` in [`packages/treetime/src/partition/marginal_sparse.rs`](../../../packages/treetime/src/partition/marginal_sparse.rs)

## The likelihood surface

The tree likelihood `L(t)` is a function of all branch lengths simultaneously. For tree refinement, we care about two aspects of this surface:

**Per-edge slices.** Fixing all branch lengths except one edge gives a 1D function of that edge's branch length. This function is what the per-edge optimizer maximizes ([Chapter 5](5-branch-length-optimization.md)). For JC69 and F81, this 1D function is unimodal. For K2P and more complex models, it can have multiple local maxima (Dinh and Matsen 2017).

**The alternating optimization landscape.** When all branch lengths are optimized in a sweep, the resulting ancestral profiles change, which shifts the optimal branch lengths for the next sweep. This creates the oscillation that damping addresses ([Chapter 9](9-iteration-loop.md)).

## References

- Felsenstein, J. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _J. Mol. Evol._ 17:368-376. https://doi.org/10.1007/BF01734359
- Felsenstein, J. 1973. "Maximum Likelihood and Minimum-Steps Methods for Estimating Evolutionary Trees." _Syst. Zool._ 22(3):240-249. https://doi.org/10.2307/2412304
- Dinh, V. C., and F. A. Matsen IV. 2017. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." _Ann. Appl. Prob._ 27(3):1646-1677. https://doi.org/10.1214/16-AAP1240
