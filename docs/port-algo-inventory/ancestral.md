# Ancestral Reconstruction Algorithms

[Back to index](_index.md)

## Fitch Parsimony

Maximum parsimony (Fitch 1971) reconstructs ancestral character states by minimizing the total number of state changes on the tree. The method makes no assumptions about branch lengths or substitution rates, treating all state transitions as equally costly. It remains widely used for seeding ML optimization with initial ancestral assignments due to its speed, and for compression of sequence data to variable-position-only representations.

v1: [`packages/treetime/src/commands/ancestral/fitch.rs#L85-L513`](../../packages/treetime/src/commands/ancestral/fitch.rs#L85-L513).
v0: [`packages/legacy/treetime/treetime/treeanc.py#L575-L686`](../../packages/legacy/treetime/treetime/treeanc.py#L575-L686).

### Algorithm

Two-pass dynamic programming on a rooted tree:

**Backward pass** (postorder, leaves to root): For each internal node, compute the set of possible states S from its children's state sets. If the children's sets intersect, S = intersection (no change needed). If they do not, S = union (at least one change occurred on the branches leading to this node). Each union operation increments the parsimony score by one.

**Forward pass** (preorder, root to leaves): Assign definite states top-down. At the root, pick one state from S_root. For each descendant, assign the parent's state if it appears in the child's state set; otherwise pick from the child's set.

### v1 extensions

- **Sparse representation**: only variable positions (positions that differ from the reference) are stored and processed. Invariant positions carry no phylogenetic signal for parsimony and can be skipped.
- **Indel handling**: insertions and deletions are tracked alongside substitutions, with majority rule for gap vs non-gap resolution at internal nodes.
- **`BitSet128` state sets**: character state sets are represented as 128-bit bitmasks, enabling O(1) intersection and union via hardware AND/OR instructions.
- **Parallel BFS traversal**: nodes at the same tree depth are processed in parallel using Rayon.

### v0 differences

v1 uses deterministic `get_one()` (`#get_one`) for root state selection when the root state set has multiple elements. v0 uses random selection (`np.random.choice`). See [intentional change](../port-intentional-changes/ancestral-fitch-deterministic-root-state.md).

### Key functions

- `fitch_backward()` (`#fitch_backward`) and `run_fitch_backward()` (`#run_fitch_backward`): postorder pass computing state sets and parsimony score
- `fitch_forward()` (`#fitch_forward`) and `run_fitch_forward()` (`#run_fitch_forward`): preorder pass resolving ambiguities

### References

- Fitch (1971). "Toward Defining the Course of Evolution: Minimum Change for a Specific Tree Topology." Systematic Zoology, 20(4):406-416. doi:10.2307/2412116
- Farris (1970). "Methods for Computing Wagner Trees." Systematic Zoology, 19(1):83-92. (Independent development of parsimony methods.)
- Felsenstein (1978). "Cases in which Parsimony or Compatibility Methods Will be Positively Misleading." Systematic Zoology, 27(4):401-410. (Long branch attraction: parsimony is statistically inconsistent when two unrelated long branches are present.)

---

## Marginal ML

Maximum likelihood ancestral reconstruction via the Felsenstein pruning algorithm (Felsenstein 1981), equivalent to the sum-product algorithm (belief propagation) on a tree-structured factor graph (Pearl 1988). Each site is treated independently: the total likelihood is a product over sites.

The algorithm computes partial likelihoods at each node - the probability of observing the data in the node's subtree given each possible ancestral state. For an internal node k with children i and j:

```
w_k(X) = [sum_Y P(X->Y|t_i) * w_i(Y)] * [sum_Z P(X->Z|t_j) * w_j(Z)]
```

where P(X->Y|t) = exp(Q*t) is the transition probability matrix from the GTR substitution model. At leaves, w is 1 for the observed state and 0 elsewhere (or spread across states for ambiguous characters). At the root, the site likelihood is `P(D_s|T) = sum_X pi_X * w_root(X)`.

The backward pass (leaf-to-root) computes partial likelihoods. The forward pass (root-to-leaf) computes "outgroup messages" via cavity/division: each node receives the information from the rest of the tree excluding its own subtree, and combines it with the backward message to produce the marginal posterior.

### v1 implementations

**Dense** (all positions): [`packages/treetime/src/representation/partition/marginal_dense.rs#L87-L281`](../../packages/treetime/src/representation/partition/marginal_dense.rs#L87-L281). Stores full probability vectors at every alignment position. Used when the full profile is needed (e.g., GTR inference from data).

**Sparse** (variable positions only): [`packages/treetime/src/representation/partition/marginal_passes.rs#L16-L250`](../../packages/treetime/src/representation/partition/marginal_passes.rs#L16-L250). Stores profiles only at positions that vary from the Fitch reference. Much faster for conserved alignments where >90% of positions are invariant.

v0: [`packages/legacy/treetime/treetime/treeanc.py#L762-L927`](../../packages/legacy/treetime/treetime/treeanc.py#L762-L927).

### v0 differences

v1 backward pass uses log-space arithmetic with logsumexp normalization (dense `normalize_from_log()`, sparse `logsumexp_normalize()` in `combine_messages()`). v1 forward pass uses plain probability space (division). v0 uses neg-log space throughout. v1 dense uses deterministic `argmax_first()` (`#argmax_first`) for sequence extraction (leftmost state wins ties); v0 uses `np.argmax()` which has undefined tie-breaking.

### Key functions

- `process_node_backward()` (`#process_node_backward`): computes partial likelihoods via GTR matrix multiplication
- `process_node_forward()` (`#process_node_forward`): computes outgroup messages via cavity/division
- `combine_messages()` (`#combine_messages`): combines child messages in log-space via logsumexp normalization
- `propagate_raw()` (`#propagate_raw`): GTR matrix-vector product `P(t) * profile`

### Complexity

O(n _ k^2 _ L) total for n nodes, k alphabet states (4 for nucleotides, 20 for amino acids), and L alignment positions. Brute force summation over all possible ancestral assignments would be O(k^n \* L), exponential in tree size.

### References

- Felsenstein (1981). "Evolutionary trees from DNA sequences: a maximum likelihood approach." J Mol Evol, 17(6):368-376. doi:10.1007/BF01734359
- Pearl (1988). "Probabilistic Reasoning in Intelligent Systems." Morgan Kaufmann. (Sum-product / belief propagation on graphical models.)
- Kschischang, Frey & Loeliger (2001). "Factor Graphs and the Sum-Product Algorithm." IEEE Trans Inform Theory, 47(2):498-519. (Formal connection between Felsenstein pruning and sum-product message passing.)

---

## Joint ML (Unimplemented)

Joint maximum likelihood reconstruction (Pupko et al. 2000) finds the single most likely assignment of ancestral states across all nodes simultaneously, rather than marginalizing over alternatives at each node independently. Uses traceback pointers (argmax) instead of marginalization (sum), analogous to the Viterbi algorithm for HMMs vs the forward-backward algorithm.

v1: `unimplemented!()` at [`packages/treetime/src/commands/ancestral/run.rs#L194`](../../packages/treetime/src/commands/ancestral/run.rs#L194). Intentionally removed - see [intentional change](../port-intentional-changes/ancestral-joint-reconstruction-removed.md).
v0: [`packages/legacy/treetime/treetime/treeanc.py#L934-L1080`](../../packages/legacy/treetime/treetime/treeanc.py#L934-L1080).

See [unimplemented](unimplemented.md#joint-ml) for full v0 algorithm details.

Reference: Pupko, Pe'er, Shamir & Graur (2000). "A fast algorithm for joint reconstruction of ancestral amino acid sequences." Mol Biol Evol, 17(6):890-896. doi:10.1093/oxfordjournals.molbev.a026369

---

## File Index

| File                                                                                                                                             | Algorithms                                                                       |
| ------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------- |
| [`packages/treetime/src/commands/ancestral/fitch.rs`](../../packages/treetime/src/commands/ancestral/fitch.rs)                                   | Fitch parsimony (backward, forward, cleanup)                                     |
| [`packages/treetime/src/commands/ancestral/marginal.rs`](../../packages/treetime/src/commands/ancestral/marginal.rs)                             | Marginal ML orchestration                                                        |
| [`packages/treetime/src/commands/ancestral/run.rs`](../../packages/treetime/src/commands/ancestral/run.rs)                                       | Ancestral command entry point, method dispatch                                   |
| [`packages/treetime/src/representation/partition/marginal_dense.rs`](../../packages/treetime/src/representation/partition/marginal_dense.rs)     | Dense marginal (Felsenstein pruning)                                             |
| [`packages/treetime/src/representation/partition/marginal_sparse.rs`](../../packages/treetime/src/representation/partition/marginal_sparse.rs)   | Sparse marginal                                                                  |
| [`packages/treetime/src/representation/partition/marginal_passes.rs`](../../packages/treetime/src/representation/partition/marginal_passes.rs)   | Sparse message passing                                                           |
| [`packages/treetime/src/representation/partition/marginal_helpers.rs`](../../packages/treetime/src/representation/partition/marginal_helpers.rs) | `combine_messages()` (`#combine_messages`), `propagate_raw()` (`#propagate_raw`) |
