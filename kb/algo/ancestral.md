# Ancestral Reconstruction Algorithms

[Back to index](README.md)

## Fitch Parsimony

Maximum parsimony (<a id="cite-1"></a>[Fitch 1971](https://doi.org/10.2307/2412116) [[1](#ref-1)]) reconstructs ancestral character states by minimizing the total number of state changes on the tree. The method makes no assumptions about branch lengths or substitution rates, treating all state transitions as equally costly. It remains widely used for seeding ML optimization with initial ancestral assignments due to its speed, and for compression of sequence data to variable-position-only representations.

v1: [`packages/treetime/src/ancestral/fitch.rs#L85-L513`](../../packages/treetime/src/ancestral/fitch.rs#L85-L513).
v0: [`packages/legacy/treetime/treetime/treeanc.py#L575-L686`](../../packages/legacy/treetime/treetime/treeanc.py#L575-L686).

### Algorithm

Two-pass dynamic programming on a rooted tree:

**Backward pass** (postorder, leaves to root): For each internal node, compute the set of possible states S from its children's state sets. If the children's sets intersect, S = intersection (no change needed). If they do not, S = union (at least one change occurred on the branches leading to this node). Each union operation increments the parsimony score by one.

**Forward pass** (preorder, root to leaves): Assign definite states top-down. At the root, pick one state from S_root. For each descendant, assign the parent's state if it appears in the child's state set; otherwise pick from the child's set.

### v1 extensions

- Sparse representation: only variable positions (positions that differ from the reference) are stored and processed. Invariant positions carry no phylogenetic signal for parsimony and can be skipped.
- Indel handling: insertions and deletions are tracked alongside substitutions, with majority rule for gap vs non-gap resolution at internal nodes.
- `BitSet128` state sets: character state sets are represented as 128-bit bitmasks, enabling O(1) intersection and union via hardware AND/OR instructions.
- Parallel BFS traversal: nodes at the same tree depth are processed in parallel using Rayon.

### v0 differences

v1 uses deterministic `get_one()` (`#get_one`) for root state selection when the root state set has multiple elements. v0 uses random selection (`np.random.choice`). See [intentional change](../decisions/ancestral-fitch-deterministic-root-state.md).

### Key functions

- `fitch_backward()` (`#fitch_backward`) and `run_fitch_backward()` (`#run_fitch_backward`): postorder pass computing state sets and parsimony score
- `fitch_forward()` (`#fitch_forward`) and `run_fitch_forward()` (`#run_fitch_forward`): preorder pass resolving ambiguities

### References

- <a id="ref-1"></a>Fitch, Walter M. 1971. "Toward Defining the Course of Evolution: Minimum Change for a Specific Tree Topology." _Systematic Zoology_ 20(4):406-416. https://doi.org/10.2307/2412116 [↩](#cite-1)
- <a id="ref-2"></a>Farris, James S. 1970. "Methods for Computing Wagner Trees." _Systematic Zoology_ 19(1):83-92. https://doi.org/10.2307/2412028
- <a id="ref-3"></a>Felsenstein, Joseph. 1978. "Cases in which Parsimony or Compatibility Methods Will be Positively Misleading." _Systematic Zoology_ 27(4):401-410. https://doi.org/10.2307/2412923

---

## Marginal ML

Maximum likelihood ancestral reconstruction via the Felsenstein pruning algorithm (<a id="cite-4"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[4](#ref-4)]), equivalent to the sum-product algorithm (belief propagation) on a tree-structured factor graph (<a id="cite-5"></a>[Pearl 1988](https://doi.org/10.1016/B978-0-08-051489-5.50001-5) [[5](#ref-5)]). Each site is treated independently: the total likelihood is a product over sites.

The algorithm computes partial likelihoods at each node - the probability of observing the data in the node's subtree given each possible ancestral state. For an internal node k with children i and j:

```
w_k(X) = [sum_Y P(X->Y|t_i) * w_i(Y)] * [sum_Z P(X->Z|t_j) * w_j(Z)]
```

where P(X->Y|t) = exp(Q*t) is the transition probability matrix from the GTR substitution model. At leaves, w is 1 for the observed state and 0 elsewhere (or spread across states for ambiguous characters). At the root, the site likelihood is `P(D_s|T) = sum_X pi_X * w_root(X)`.

The backward pass (leaf-to-root) computes partial likelihoods. The forward pass (root-to-leaf) computes "outgroup messages" via cavity/division: each node receives the information from the rest of the tree excluding its own subtree, and combines it with the backward message to produce the marginal posterior.

### v1 implementations

**Dense** (all positions): [`packages/treetime/src/partition/marginal_dense.rs#L87-L281`](../../packages/treetime/src/partition/marginal_dense.rs#L87-L281). Stores full probability vectors at every alignment position. Used when the full profile is needed (e.g., GTR inference from data).

**Sparse** (variable positions only): [`packages/treetime/src/partition/marginal_passes.rs#L16-L250`](../../packages/treetime/src/partition/marginal_passes.rs#L16-L250). Stores profiles only at positions that vary from the Fitch reference. Much faster for conserved alignments where >90% of positions are invariant.

v0: [`packages/legacy/treetime/treetime/treeanc.py#L762-L927`](../../packages/legacy/treetime/treetime/treeanc.py#L762-L927).

### v0 differences

v1 backward pass uses log-space arithmetic with logsumexp normalization (dense `normalize_from_log()`, sparse `softmax_with_log_norm()` in `combine_messages()`). v1 forward pass uses plain probability space (division). v0 uses neg-log space throughout. v1 dense uses deterministic `argmax_first()` (`#argmax_first`) for sequence extraction (leftmost state wins ties); v0 uses `np.argmax()` which has undefined tie-breaking.

### Key functions

- `process_node_backward()` (`#process_node_backward`): computes partial likelihoods via GTR matrix multiplication
- `process_node_forward()` (`#process_node_forward`): computes outgroup messages via cavity/division
- `combine_messages()` (`#combine_messages`): combines child messages in log-space via logsumexp normalization
- `propagate_raw()` (`#propagate_raw`): GTR matrix-vector product `P(t) * profile`

### Complexity

O(n _ k^2 _ L) total for n nodes, k alphabet states (4 for nucleotides, 20 for amino acids), and L alignment positions. Brute force summation over all possible ancestral assignments would be O(k^n \* L), exponential in tree size.

### References

- <a id="ref-4"></a>Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-4)
- <a id="ref-5"></a>Pearl, Judea. 1988. _Probabilistic Reasoning in Intelligent Systems: Networks of Plausible Inference._ Morgan Kaufmann. ISBN 978-0-934613-73-2. [↩](#cite-5)
- <a id="ref-6"></a>Kschischang, Frank R., Brendan J. Frey, and Hans-Andrea Loeliger. 2001. "Factor Graphs and the Sum-Product Algorithm." _IEEE Transactions on Information Theory_ 47(2):498-519. https://doi.org/10.1109/18.910572

---

## Joint ML (Unimplemented)

Joint maximum likelihood reconstruction (<a id="cite-7"></a>[Pupko et al. 2000](https://doi.org/10.1093/oxfordjournals.molbev.a026369) [[7](#ref-7)]) finds the single most likely assignment of ancestral states across all nodes simultaneously, rather than marginalizing over alternatives at each node independently. Uses traceback pointers (argmax) instead of marginalization (sum), analogous to the Viterbi algorithm for HMMs vs the forward-backward algorithm.

v1: `unimplemented!()` at [`packages/treetime/src/commands/ancestral/run.rs#L199`](../../packages/treetime/src/commands/ancestral/run.rs#L199). Intentionally removed - see [intentional change](../decisions/ancestral-joint-reconstruction-removed.md).
v0: [`packages/legacy/treetime/treetime/treeanc.py#L934-L1080`](../../packages/legacy/treetime/treetime/treeanc.py#L934-L1080).

See [unimplemented](unimplemented.md#joint-ml) for full v0 algorithm details.

### References

- <a id="ref-7"></a>Pupko, Tal, Itsik Pe'er, Ron Shamir, and Dan Graur. 2000. "A Fast Algorithm for Joint Reconstruction of Ancestral Amino Acid Sequences." _Molecular Biology and Evolution_ 17(6):890-896. https://doi.org/10.1093/oxfordjournals.molbev.a026369 [↩](#cite-7)

---

## Branch Mutation Annotation

After ancestral reconstruction (Fitch or marginal), branch mutations are extracted from partition data and attached to tree nodes for Newick/Nexus output. The annotation pipeline is shared across all commands that output annotated trees (ancestral, timetree, optimize).

v1: [`packages/treetime/src/payload/ancestral.rs#L152-L186`](../../packages/treetime/src/payload/ancestral.rs#L152-L186) (`annotate_branch_mutations()`).
v0: annotation is inline in `treeanc.py` tree-writing methods.

### Algorithm

For each edge in the tree, `annotate_branch_mutations()` calls `PartitionBranchOps::edge_subs()` on every partition, collects all substitutions, sorts them by position, and writes the formatted comma-separated string (e.g. `"A55G,T93C"`, 1-based positions) into the child node's `mutations` field. This field is emitted as a `mutations="..."` NHX comment in Newick/Nexus output.

Both dense and sparse `edge_subs()` implementations apply the same filtering: only canonical nucleotide substitutions where parent and child states differ are reported. Gaps, unknowns, and ambiguity codes are excluded.

### Dense `edge_subs()`

[`packages/treetime/src/partition/marginal_dense.rs#L81-L112`](../../packages/treetime/src/partition/marginal_dense.rs#L81-L112)

Iterates every alignment position (0..L where L = number of rows in the profile matrix). At each position, takes the MAP state (`argmax_first()` of the posterior profile) for both parent and child. Skips positions that fall within either endpoint's original gap ranges (`DenseSeqInfo.gaps`), because gap positions receive uniform profiles under `treat_gap_as_unknown` and their argmax would return an arbitrary canonical state.

### Sparse `edge_subs()`

[`packages/treetime/src/partition/marginal_sparse.rs#L233-L239`](../../packages/treetime/src/partition/marginal_sparse.rs#L233-L239)

Returns MAP-derived substitutions stored in `subs_ml`. Requires marginal inference to have run (errors if `subs_ml` is `None`). The ML subs are computed during the marginal forward pass by `compute_ml_subs_for_edge()` in [`marginal_passes.rs`](../../packages/treetime/src/partition/marginal_passes.rs), which compares parent and child MAP states at candidate positions (union of Fitch subs and variable sites). ML subs are cleared automatically by any fitch-sub mutation and by `clear_ml_subs()` during reroot.

The candidate set is complete: any position where parent and child could differ must appear as a variable site on at least one endpoint or as a Fitch substitution on the edge.

### Dense vs sparse equivalence

Both implementations produce the same mutation set for the same reconstruction. Dense scans all L positions but most comparisons are equal (no-op). Sparse scans only the variable-site union, which for conserved alignments (>90% invariant) is a small fraction of L.

### Call sites

- [`packages/treetime/src/commands/ancestral/run.rs#L142`](../../packages/treetime/src/commands/ancestral/run.rs#L142) and [`#L192`](../../packages/treetime/src/commands/ancestral/run.rs#L192): ancestral command (both dense and sparse paths)
- [`packages/treetime/src/commands/timetree/run.rs#L490`](../../packages/treetime/src/commands/timetree/run.rs#L490): timetree command
- [`packages/treetime/src/commands/optimize/run.rs#L225`](../../packages/treetime/src/commands/optimize/run.rs#L225): optimize command

### Key functions

- `annotate_branch_mutations()` (`#annotate_branch_mutations`): generic over graph payload, iterates edges and partitions, writes formatted mutation string to child node
- `PartitionBranchOps::edge_subs()` (`#edge_subs`): trait method implemented by both `PartitionMarginalDense` and `PartitionMarginalSparse`
- `compute_ml_subs_for_edge()` (`#compute_ml_subs_for_edge`): sparse-only, computes MAP-derived subs from finalized parent and child profiles during the forward pass, stores result in `subs_ml`
- `reconstruct_map_sequence()` (`#reconstruct_map_sequence`): sparse-only, rebuilds node sequence from parent + edge mutations + MAP variable-site states during the forward pass
- `HasBranchMutations::set_branch_mutations()` (`#set_branch_mutations`): trait implemented by `NodeAncestral` and `NodeTimetree` (via delegation to inner `NodeAncestral`)

---

## File Index

| File                                                                                                                                             | Algorithms                                                                       |
| ------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------- |
| [`packages/treetime/src/ancestral/fitch.rs`](../../packages/treetime/src/ancestral/fitch.rs)                                   | Fitch parsimony (backward, forward, cleanup)                                     |
| [`packages/treetime/src/ancestral/marginal.rs`](../../packages/treetime/src/ancestral/marginal.rs)                             | Marginal ML orchestration                                                        |
| [`packages/treetime/src/commands/ancestral/run.rs`](../../packages/treetime/src/commands/ancestral/run.rs)                                       | Ancestral command entry point, method dispatch                                   |
| [`packages/treetime/src/partition/marginal_dense.rs`](../../packages/treetime/src/partition/marginal_dense.rs)     | Dense marginal (Felsenstein pruning)                                             |
| [`packages/treetime/src/partition/marginal_sparse.rs`](../../packages/treetime/src/partition/marginal_sparse.rs)   | Sparse marginal                                                                  |
| [`packages/treetime/src/partition/marginal_passes.rs`](../../packages/treetime/src/partition/marginal_passes.rs)   | Sparse message passing                                                           |
| [`packages/treetime/src/partition/marginal_helpers.rs`](../../packages/treetime/src/partition/marginal_helpers.rs) | `combine_messages()` (`#combine_messages`), `propagate_raw()` (`#propagate_raw`) |
| [`packages/treetime/src/payload/ancestral.rs`](../../packages/treetime/src/payload/ancestral.rs)                   | Branch mutation annotation (`annotate_branch_mutations()`)                       |
| [`packages/treetime/src/partition/traits.rs`](../../packages/treetime/src/partition/traits.rs)                     | `PartitionBranchOps` trait (`edge_subs()`)                                       |
