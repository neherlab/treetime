# Chapter 4: Ancestral sequence reconstruction

[Back to index](_index.md) | Previous: [Chapter 3: Tree likelihood](3-tree-likelihood.md) | Next: [Chapter 5: Branch length optimization](5-branch-length-optimization.md)

## What is ancestral reconstruction?

Ancestral sequence reconstruction infers the nucleotide sequences at internal nodes of the phylogenetic tree. These ancestral sequences are not observed -- only the leaf sequences (sampled organisms) are known. The reconstruction fills in the missing data, assigning a nucleotide (or a probability distribution over nucleotides) to each internal node at each alignment position.

Ancestral reconstruction serves two purposes in tree refinement:

1. **Branch length optimization** ([Chapter 5](5-branch-length-optimization.md)): the per-edge optimizer needs to know what states are at both endpoints of each branch. The ancestral reconstruction provides the internal-node endpoints.

2. **Topology cleanup**: the `prune_short_branches()` criterion in v0 evaluates the probability of the parent-child sequence pair at zero distance ([Chapter 6](6-zero-length-branches.md)). The shared-mutation merging algorithm compares substitution sets on sibling branches ([Chapter 7](7-polytomy-resolution.md)). Both require knowing which mutations occurred on which branches, which requires ancestral reconstruction.

Two families of methods exist: parsimony (fast, approximate) and maximum likelihood (slower, exact).

## Parsimony: Fitch's algorithm

Fitch (1971) introduced the standard parsimony algorithm for ancestral state reconstruction. The algorithm minimizes the total number of state changes on the tree, treating all substitutions as equally costly. It makes no assumptions about branch lengths or substitution rates.

### The algorithm

Two-pass dynamic programming on a rooted tree:

**Backward pass** (postorder, leaves to root): At each internal node, compute the **state set** S -- the set of states that minimize the parsimony cost in the subtree below. If the children's sets intersect, `S = intersection` (no change needed at this node). If they do not intersect, `S = union` (at least one change occurred), and the parsimony score increments by one.

```
Example (position i):

  Leaf A: {G}     Leaf B: {G}     Leaf C: {T}     Leaf D: {T}

  Node AB: {G} intersect {G} = {G}    (no change)
  Node CD: {T} intersect {T} = {T}    (no change)
  Root: {G} intersect {T} = empty, so {G, T} union  (1 change)
```

**Forward pass** (preorder, root to leaves): Assign definite states top-down. At the root, pick one state from the root's set. For each descendant, use the parent's assigned state if it appears in the child's set; otherwise pick from the child's set.

```
Forward pass (continuing example, root picks G):

  Root: G
  Node AB: G in {G} -> G
  Node CD: G not in {T} -> T    (change on branch Root-CD)
  Leaf A: G    Leaf B: G    Leaf C: T    Leaf D: T
```

### Complexity

O(n \* s) per site, where n is the number of nodes and s is the number of states. For nucleotides (s=4) with bit-set representation, intersection and union are single AND/OR instructions. TreeTime v1 uses `BitSet128` for O(1) set operations.

### Role in TreeTime

Fitch reconstruction is used for:

- **Compression**: identifying variable positions (positions where the state set at the root has more than one element). Invariant positions are skipped in subsequent ML computation.
- **Initial ancestral assignment**: seeding the ML optimization with a parsimony-based starting point.
- **Mutation mapping**: determining which branches carry substitutions, used by `--prune-empty` and `--merge-shared-mutations`.

v1 code: [`packages/treetime/src/commands/ancestral/fitch.rs`](../../../packages/treetime/src/commands/ancestral/fitch.rs). The `compress_sequences()` function runs Fitch reconstruction and stores variable-site data in the sparse partition structure.

v0 code: [`packages/legacy/treetime/treetime/treeanc.py#L575-L686`](../../../packages/legacy/treetime/treetime/treeanc.py#L575-L686). The `_fitch_anc()` method.

### Limitations

Parsimony does not account for branch lengths or substitution rates. On trees with long branches, it can produce incorrect reconstructions due to **long-branch attraction**: convergent substitutions on long branches are mistaken for shared ancestry. For tree refinement, parsimony is used only as an initial approximation; the ML methods below provide the definitive reconstruction.

## Maximum likelihood: marginal reconstruction

Marginal reconstruction computes the posterior probability of each state at each internal node independently, given all the data in the tree. At each node u and position i, the result is a probability vector:

```
P(state_u = a | D) = [P(A), P(C), P(G), P(T)]
```

This is the product of two messages from the pruning algorithm ([Chapter 3](3-tree-likelihood.md)):

```
P(state_u = a | D) proportional to msg_to_parent(u)[a] * msg_to_child(u)[a]
```

where `msg_to_parent` carries information from the subtree below u, and `msg_to_child` carries information from the rest of the tree above u.

### The backward pass (subtree messages)

Same as the pruning algorithm's likelihood computation ([Chapter 3](3-tree-likelihood.md)): traverse from leaves to root, computing conditional likelihood vectors at each node. The result at each node is the probability of the subtree data given each possible state at that node.

### The forward pass (outgroup messages)

Traverse from root to leaves. At each node u with parent p, the outgroup message combines the parent's outgroup message with the sibling subtree messages:

```
msg_to_child(u)[a] = sum_b [msg_to_child(p)[b] * prod_{siblings s of u} msg_to_parent(s)[b]] * P(t_u)_{ba}
```

This tells node u: "given everything the rest of the tree says about what state you should have, here is the probability of each state."

### The MAP state

The **maximum a posteriori (MAP)** state at each position is the argmax of the posterior:

```
MAP_u(i) = argmax_a P(state_u = a | D)
```

The MAP states define the reconstructed ancestral sequence. Branch substitutions are differences between the MAP states at parent and child nodes. This is how v1's `edge_subs()` function determines which mutations occurred on each branch.

v1 code:

- Backward pass: `process_node_backward()` in [`packages/treetime/src/representation/partition/marginal_passes.rs`](../../../packages/treetime/src/representation/partition/marginal_passes.rs)
- Forward pass: `process_node_forward()` in the same file
- Marginal orchestration: `update_marginal()` in [`packages/treetime/src/commands/ancestral/marginal.rs`](../../../packages/treetime/src/commands/ancestral/marginal.rs)
- Dense edge substitutions: `edge_subs()` on `PartitionMarginalDense` compares MAP states at parent and child node posteriors
- Sparse edge substitutions: `edge_subs()` on `PartitionMarginalSparse` returns marginal-reconstructed substitutions

v0 code: `_ml_anc_marginal()` in [`packages/legacy/treetime/treetime/treeanc.py`](../../../packages/legacy/treetime/treetime/treeanc.py).

### Marginal vs joint reconstruction

**Marginal** reconstruction optimizes each position independently. The reconstruction at position i does not depend on the reconstruction at position j. This is correct under the standard assumption that sites evolve independently.

**Joint** reconstruction finds the single most likely assignment of states to all internal nodes simultaneously. Joint reconstruction was implemented in v0 using the Viterbi-like algorithm of Pupko et al. (2000) but has been removed in v1 as an intentional simplification. Joint reconstruction can produce different ancestral sequences than marginal reconstruction at positions where the posterior is multimodal (multiple states have similar probabilities).

For tree refinement, marginal reconstruction is standard because:

- It provides probability distributions (not just point estimates), enabling soft comparisons
- It is the natural E-step of the EM algorithm ([Chapter 9](9-iteration-loop.md))
- It avoids the computational cost of joint optimization across all nodes

## Interaction with branch length optimization

Ancestral reconstruction and branch length optimization are interdependent:

1. **Reconstruction depends on branch lengths.** The transition probability matrices `P(t) = exp(Q*t)` in the pruning algorithm use the current branch lengths. Different branch lengths produce different posterior profiles at internal nodes.

2. **Optimization depends on reconstruction.** The per-edge optimizer ([Chapter 5](5-branch-length-optimization.md)) uses the `msg_to_parent` and `msg_to_child` messages to compute the likelihood and its derivatives. These messages come from the reconstruction.

This circular dependence is why tree refinement requires an iterative loop ([Chapter 9](9-iteration-loop.md)): reconstruct given current branch lengths, optimize branch lengths given current reconstruction, repeat.

## References

- Fitch, W. M. 1971. "Toward Defining the Course of Evolution: Minimum Change for a Specific Tree Topology." _Syst. Zool._ 20(4):406-416. https://doi.org/10.2307/2412116
- Felsenstein, J. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _J. Mol. Evol._ 17:368-376. https://doi.org/10.1007/BF01734359
- Pupko, T., I. Pe'er, R. Shamir, and D. Graur. 2000. "A Fast Algorithm for Joint Reconstruction of Ancestral Amino Acid Sequences." _Mol. Biol. Evol._ 17(6):890-896. https://doi.org/10.1093/oxfordjournals.molbev.a026369
- Yang, Z., S. Kumar, and M. Nei. 1995. "A New Method of Inference of Ancestral Nucleotide and Amino Acid Sequences." _Genetics_ 141(4):1641-1650. https://doi.org/10.1093/genetics/141.4.1641
