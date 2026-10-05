# Ancestral Reconstruction Algorithms

## Fitch Parsimony

Maximum parsimony (<a id="cite-1"></a>[Fitch 1971](https://doi.org/10.2307/2412116) [[1](#ref-1)]) reconstructs ancestral character states by minimizing the total number of state changes on the tree. The method makes no assumptions about branch lengths or substitution rates, treating all state transitions as equally costly. It remains widely used for seeding ML optimization with initial ancestral assignments due to its speed, and for compression of sequence data to variable-position-only representations.

v1: `fn fitch_backward()`, `fn fitch_forward()`, and `fn fitch_cleanup()` in [`packages/treetime/src/partition/fitch/passes.rs#L91-L277`](../../packages/treetime/src/partition/fitch/passes.rs#L91-L277), driven per partition by [`packages/treetime/src/ancestral/fitch.rs`](../../packages/treetime/src/ancestral/fitch.rs).
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

- `fitch_backward()` (`#fitch_backward`) and `run_fitch_backward_indexed()` (`#run_fitch_backward_indexed`): postorder pass computing state sets and parsimony score
- `fitch_forward()` (`#fitch_forward`) and `run_fitch_forward_indexed()` (`#run_fitch_forward_indexed`): preorder pass resolving ambiguities

### References

- <a id="ref-1"></a>Fitch, Walter M. 1971. "Toward Defining the Course of Evolution: Minimum Change for a Specific Tree Topology." _Systematic Zoology_ 20(4):406-416. https://doi.org/10.2307/2412116 [↩](#cite-1)
- <a id="ref-2"></a>Farris, James S. 1970. "Methods for Computing Wagner Trees." _Systematic Zoology_ 19(1):83-92. https://doi.org/10.2307/2412028
- <a id="ref-3"></a>Felsenstein, Joseph. 1978. "Cases in which Parsimony or Compatibility Methods Will be Positively Misleading." _Systematic Zoology_ 27(4):401-410. https://doi.org/10.2307/2412923

---

## Marginal ML

Maximum likelihood ancestral reconstruction via the Felsenstein pruning algorithm (<a id="cite-4"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[4](#ref-4)]), equivalent to the sum-product algorithm (belief propagation) on a tree-structured factor graph (<a id="cite-5"></a>[Pearl 1988](https://doi.org/10.1016/C2009-0-27609-4) [[5](#ref-5)]). Each site is treated independently: the total likelihood is a product over sites.

The algorithm computes partial likelihoods at each node - the probability of observing the data in the node's subtree given each possible ancestral state. For an internal node k with children i and j:

```
w_k(X) = [sum_Y P(X->Y|t_i) * w_i(Y)] * [sum_Z P(X->Z|t_j) * w_j(Z)]
```

where P(X->Y|t) = exp(Q*t) is the transition probability matrix from the GTR substitution model. At leaves, w is 1 for the observed state and 0 elsewhere (or spread across states for ambiguous characters). At the root, the site likelihood is `P(D_s|T) = sum_X pi_X * w_root(X)`.

The backward pass (leaf-to-root) computes partial likelihoods. The forward pass (root-to-leaf) computes "outgroup messages" via cavity/division: each node receives the information from the rest of the tree excluding its own subtree, and combines it with the backward message to produce the marginal posterior.

Both passes use [`packages/treetime-graph/src/pass.rs`](../../packages/treetime-graph/src/pass.rs) to move partition maps into pass-scoped topology-indexed storage. Each node owns a disjoint mutable slot. Backward work reads only completed descendant payloads; forward work reads only completed ancestor payloads. Parent-edge data has a stable slot throughout the pass. A node runs as soon as its dependencies complete rather than at a per-level barrier, preserving the sequential dependency order while allowing independent nodes to run concurrently.

### v1 implementations

**Dense** (all positions): `fn indexed_backward()` and `fn indexed_forward()` in [`packages/treetime/src/partition/marginal/shared/pass.rs`](../../packages/treetime/src/partition/marginal/shared/pass.rs) (shared with discrete partitions), partition in [`packages/treetime/src/partition/marginal/dense/partition.rs`](../../packages/treetime/src/partition/marginal/dense/partition.rs). Stores full probability vectors at every alignment position. Used when the full profile is needed (e.g., GTR inference from data).

**Sparse** (variable positions only): [`packages/treetime/src/partition/marginal/sparse/`](../../packages/treetime/src/partition/marginal/sparse/) (`backward.rs`, `forward.rs`, `message.rs`). Stores explicit profiles only at positions that vary near each node; see [Sparse marginal design](#sparse-marginal-design). Much faster for conserved alignments where >90% of positions are invariant.

v0: [`packages/legacy/treetime/treetime/treeanc.py#L762-L927`](../../packages/legacy/treetime/treetime/treeanc.py#L762-L927).

### Sparse node data model

Each sparse node has observations (`struct SparseNodeObs`) and a pass state (`struct SparseNodeState`), both in [`packages/treetime/src/partition/storage/sparse.rs`](../../packages/treetime/src/partition/storage/sparse.rs):

Observations: `composition` (character counts), `gaps`/`unknown`/`non_char` (position ranges where the node has no informative data), and the Fitch state. MAP reconstruction replaces the node sequence and recomputes `composition` atomically, so later sparse passes use counts from the current reconstruction rather than the earlier Fitch assignment.

Fitch state: `fitch.variable` (ambiguous positions from parsimony, as state sets), `fitch.chosen_state` (resolved state per variable position from Fitch forward pass). After the Fitch passes, `fn fitch_cleanup()` [`packages/treetime/src/partition/fitch/passes.rs#L270-L277`](../../packages/treetime/src/partition/fitch/passes.rs#L270-L277) clears `fitch.variable` for internal nodes, so only leaves keep it.

Pass state: `sequence` (stored sequence) and `profile` (`struct SparseSeqDistribution`). A profile has two parts. `variable` maps each explicit position to a `VarPos { dis, state }`: `dis` is the probability vector of the position and `state` its Fitch state. `fixed` maps each canonical state to one shared vector, used for every position not listed in `variable`, and `fixed_counts` counts those positions per state. A position that is certainly `A` in every message needs no vector of its own; it is one count in `fixed_counts['A']` and uses `fixed['A']`.

`composition` flows into `fixed_counts` in `msg_to_child` (`fn compute_msg_to_child()`), which provides the per-character multiplicity weights for branch-length optimization ([`packages/treetime/src/partition/optimize/sparse.rs`](../../packages/treetime/src/partition/optimize/sparse.rs)). `non_char` is used to mask positions from state propagation in both passes and to exclude non-evolving positions from `edge_effective_length`.

Non-char (N, gap) differences between parent and child are tracked through `non_char` ranges, not as Fitch substitutions. Applying a child edge's fitch_subs and indels to the parent's composition does not reproduce the child's composition when non-char positions differ. This is by design: Fitch parsimony assigns concrete states at ambiguous positions, and N/gap masking is a separate layer.

### Sparse marginal design

Branches in typical TreeTime datasets are short, and on a given edge or near a given node only a small minority of sites change. Computing propagators for every site wastes most of the work on sites with negligible probability of change. Two alternatives exist. Parsimony avoids probabilistic models altogether. MAPLE (<a id="cite-8"></a>[De Maio et al. 2023](https://doi.org/10.1038/s41588-023-01368-0) [[8](#ref-8)]) stores explicit vectors only for positions that differ from a reference genome, but that set grows with the distance from the reference, so the representation stops being sparse as the tree gets deeper. v1 instead keeps a site explicit only where it varies near the current node: at Fitch substitutions on adjacent edges, at ambiguous leaf sites, and where the messages disagree or are uncertain. Every other site is represented by the shared `fixed` vector of its Fitch state. The Fitch reconstruction supplies everything the passes need: the root sequence, the substitutions on each edge, the gap and unknown ranges of each node, the ambiguous sites of each leaf with their chosen states, and the state counts of each node.

**Backward pass** (`fn process_node_backward_indexed()` [`packages/treetime/src/partition/marginal/sparse/backward.rs#L45-L168`](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L45-L168)):

- Leaf: `fixed` holds the one-hot profile of each determined state, and `fixed_counts` is the leaf composition. Each ambiguous position becomes a `VarPos` whose `dis` is the profile of its ambiguity set and whose `state` is the Fitch chosen state (or one member of the set when the chosen state is not canonical)
- Internal node, collection: the explicit positions of the node are (1) every position with a Fitch substitution on a child edge, with the substitution's reference state as the node state and its query state recorded for that child, and (2) every position explicit in a child message, with the child's state. Substitutions take precedence
- Internal node, child states: for each child and each collected position without a substitution on that child's edge, the child state is the node state, except at the child's `non_char` positions, where it is gap or unknown. The child state selects the `fixed` vector used when the position is not explicit in that child's message. This is why `combine_messages()` takes per-message states (`reference_states`): the same function serves the forward pass, where the two messages come from the parent and the child, so the name `child_states` would be wrong
- `combine_messages()` then multiplies the messages and propagates the result to the parent edge with $P(t)^T$. At the root it also multiplies in the equilibrium frequencies $\pi$

**Message combination** (`fn combine_messages()` [`packages/treetime/src/partition/marginal/sparse/message.rs#L22-L103`](../../packages/treetime/src/partition/marginal/sparse/message.rs#L22-L103)): for each collected position, it sums the log of the message vectors in log space (the explicit vector of a message, or the `fixed` vector of that message's state; non-canonical states contribute nothing) and normalizes with `softmax_with_log_norm_owned()`. A position stays explicit unless both conditions hold:

- the largest probability exceeds $1 - \epsilon$ with $\epsilon = 10^{-4}$ (`fn is_site_resolved()`)
- every message agrees with the node's Fitch state at the position

The check does not require the argmax to equal the Fitch state; both conditions together make a disagreement unlikely but not impossible. A dropped position returns to `fixed_counts`. The `fixed` vectors are combined once per state and weighted by the number of positions that remain fixed. Dropping positions is an approximation; see [kb/decisions/ancestral-marginal-sparse-parsimony-chain.md](../decisions/ancestral-marginal-sparse-parsimony-chain.md) and [kb/issues/M-ancestral-sparse-dense-internal-residues-diverge.md](../issues/M-ancestral-sparse-dense-internal-residues-diverge.md).

**Forward pass** (`fn process_node_forward_indexed()` [`packages/treetime/src/partition/marginal/sparse/forward.rs#L69-L164`](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L69-L164)):

- `msg_to_child` (`fn compute_msg_to_child()` [`packages/treetime/src/partition/marginal/sparse/forward.rs#L173-L240`](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L173-L240)) divides the parent profile by the child's `msg_from_child`, position by position and for each `fixed` vector, and normalizes. Its explicit positions are the edge substitutions and the explicit positions of the parent profile and of `msg_from_child`, without the child's `non_char` positions
- `msg_to_child` is propagated along the edge with $P(t)$ to give `msg_from_parent`
- The node profile is `combine_messages(msg_from_parent, msg_to_parent)` without a root prior. Its explicit positions are the edge substitutions (child state), the explicit positions of `msg_to_child` outside the child's `non_char` ranges, and the explicit positions of `msg_to_parent`

Both passes also handle per-site rates (`fn propagate_raw_per_site()`), the `transmission` ranges of an edge, and the log-likelihood bookkeeping of the division step.

### v0 differences

v1 backward pass uses log-space arithmetic with logsumexp normalization (dense `normalize_from_log()`, sparse `softmax_with_log_norm()` in `combine_messages()`). v1 forward pass uses plain probability space (division). v0 uses neg-log space throughout. v1 dense uses deterministic `argmax_first()` (`#argmax_first`) for sequence extraction (leftmost state wins ties); v0 uses `np.argmax()` which has undefined tie-breaking.

### Key functions

- `process_node_backward_indexed()` (`#process_node_backward_indexed`, sparse) and `indexed_node_backward()` (`#indexed_node_backward`, dense): compute partial likelihoods via GTR matrix multiplication
- `process_node_forward_indexed()` (`#process_node_forward_indexed`, sparse) and `indexed_node_forward()` (`#indexed_node_forward`, dense): compute outgroup messages via cavity/division
- `combine_messages()` (`#combine_messages`): combines child messages in log-space via logsumexp normalization
- `propagate_raw()` (`#propagate_raw`): GTR matrix-vector product `P(t) * profile`

### Complexity

O(n _ k^2 _ L) total for n nodes, k alphabet states (4 for nucleotides, 20 for amino acids), and L alignment positions. Brute force summation over all possible ancestral assignments would be O(k^n \* L), exponential in tree size.

### References

- <a id="ref-4"></a>Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-4)
- <a id="ref-5"></a>Pearl, Judea. 1988. _Probabilistic Reasoning in Intelligent Systems: Networks of Plausible Inference._ Morgan Kaufmann. ISBN 978-0-934613-73-2. https://doi.org/10.1016/C2009-0-27609-4 [↩](#cite-5)
- <a id="ref-6"></a>Kschischang, Frank R., Brendan J. Frey, and Hans-Andrea Loeliger. 2001. "Factor Graphs and the Sum-Product Algorithm." _IEEE Transactions on Information Theory_ 47(2):498-519. https://doi.org/10.1109/18.910572
- <a id="ref-8"></a>De Maio, Nicola, Prabhav Kalaghatgi, Yatish Turakhia, Russell Corbett-Detig, Bui Quang Minh, and Nick Goldman. 2023. "Maximum Likelihood Pandemic-Scale Phylogenetics." _Nature Genetics_ 55(5):746-752. https://doi.org/10.1038/s41588-023-01368-0 [↩](#cite-8)

---

## Joint ML (Unimplemented)

Joint maximum likelihood reconstruction (<a id="cite-7"></a>[Pupko et al. 2000](https://doi.org/10.1093/oxfordjournals.molbev.a026369) [[7](#ref-7)]) finds the single most likely assignment of ancestral states across all nodes simultaneously, rather than marginalizing over alternatives at each node independently. Uses traceback pointers (argmax) instead of marginalization (sum), analogous to the Viterbi algorithm for HMMs vs the forward-backward algorithm.

v1: no joint mode; `--method-anc` offers no joint method. Intentionally removed - see [intentional change](../decisions/ancestral-joint-reconstruction-removed.md).
v0: [`packages/legacy/treetime/treetime/treeanc.py#L934-L1080`](../../packages/legacy/treetime/treetime/treeanc.py#L934-L1080).

See [unimplemented](unimplemented.md#joint-ml) for full v0 algorithm details.

### References

- <a id="ref-7"></a>Pupko, Tal, Itsik Pe'er, Ron Shamir, and Dan Graur. 2000. "A Fast Algorithm for Joint Reconstruction of Ancestral Amino Acid Sequences." _Molecular Biology and Evolution_ 17(6):890-896. https://doi.org/10.1093/oxfordjournals.molbev.a026369 [↩](#cite-7)

---

## Branch Mutations

The `ancestral` and `timetree` commands write the substitutions between the emitted parent and child sequences of every edge, together with the edge's indels. The `optimize` command writes the substitutions of the marginal reconstruction engine instead ([kb/issues/M-optimize-mutations-not-derived-from-emitted-sequences.md](../issues/M-optimize-mutations-not-derived-from-emitted-sequences.md)), and the `prune` command writes the Fitch substitutions of its partition.

v1: `MarginalReconstruction::stream_sequences()` in [`packages/treetime/src/partition/marginal/reconstruction.rs`](../../packages/treetime/src/partition/marginal/reconstruction.rs), called by `AncestralPartition::stream_sequences()` in [`packages/treetime/src/ancestral/partition.rs`](../../packages/treetime/src/ancestral/partition.rs) and by `fn reconstruct_final_sequences()` in [`packages/treetime/src/timetree/pipeline.rs`](../../packages/treetime/src/timetree/pipeline.rs).

### Rule

At every position, a substitution is reported when the parent and child states differ and neither is a gap. A substitution into `N` on the edge into a leaf is derived only with `--report-ambiguous`, because no descendant can bridge it. Other substitutions to and from `N` are derived and then resolved by `UnknownMutationFilter` in [`packages/app-output/src/mutation_filter.rs`](../../packages/app-output/src/mutation_filter.rs): unless `--report-ambiguous` is set, it drops them and reports a change across a masked stretch on the edge where the residue is observed again.

### Full sequence comparison

`fn stream_sequence_mutations()` in [`packages/treetime/src/seq/mutation.rs`](../../packages/treetime/src/seq/mutation.rs) materializes every node sequence in preorder, writes it to the sequence sink, and compares each child with its parent over all positions (`fn sequence_subs()`). It serves dense reconstructions, Fitch reconstructions, and marginal reconstructions with sampled sequences (`--sample-from-profile`).

### Sparse derivation

`fn sparse_edge_mutations()` in [`packages/treetime/src/partition/marginal/sparse/mutations.rs`](../../packages/treetime/src/partition/marginal/sparse/mutations.rs) gives the same result for sparse marginal reconstructions without materializing the sequences. Each node state holds a stored sequence and a set of variable positions. The emitted sequence of an internal node differs from its stored sequence only at its variable positions; a leaf's emitted sequence differs from its stored sequence only at its ambiguous positions and, with imputation, its unknown positions. So the parent and child emitted states can differ only at:

- positions where the two stored sequences differ, found by comparing them in 64-byte blocks
- the parent's variable positions
- the child's variable positions (internal child) or its ambiguous and unknown positions (leaf)

The function evaluates both emitted states at these positions only, with the same per-position rules as the full comparison: the argmax of the variable profile for internal nodes, and `fn impute_state()` in [`packages/treetime/src/partition/marginal/sparse/reconstruct.rs`](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs) for leaves. The candidate set does not depend on the edge's Fitch substitutions, because topology changes such as polytomy resolution create edges without Fitch data.

The property test `test_prop_sparse_edge_mutations_match_full_sequence_comparison` in [`packages/treetime/src/ancestral/__tests__/test_sparse_mutations_prop.rs`](../../packages/treetime/src/ancestral/__tests__/test_sparse_mutations_prop.rs) checks equality with the full comparison, with and without imputation and with the edge Fitch substitutions removed.

---

## File Index

| File                                                                                                                                 | Algorithms                                                                       |
| ------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------- |
| [`packages/treetime/src/partition/fitch/passes.rs`](../../packages/treetime/src/partition/fitch/passes.rs)                           | Fitch parsimony (backward, forward, cleanup)                                     |
| [`packages/treetime/src/ancestral/fitch.rs`](../../packages/treetime/src/ancestral/fitch.rs)                                         | Fitch reconstruction over all partitions                                         |
| [`packages/treetime/src/ancestral/pipeline.rs`](../../packages/treetime/src/ancestral/pipeline.rs)                                   | Ancestral pipeline, method dispatch                                              |
| [`packages/app-commands/src/commands/ancestral/run.rs`](../../packages/app-commands/src/commands/ancestral/run.rs)                   | Ancestral command entry point                                                    |
| [`packages/treetime/src/partition/marginal/shared/pass.rs`](../../packages/treetime/src/partition/marginal/shared/pass.rs)           | Dense and discrete marginal passes (Felsenstein pruning)                         |
| [`packages/treetime/src/partition/marginal/sparse/backward.rs`](../../packages/treetime/src/partition/marginal/sparse/backward.rs)   | Sparse marginal backward pass                                                    |
| [`packages/treetime/src/partition/marginal/sparse/forward.rs`](../../packages/treetime/src/partition/marginal/sparse/forward.rs)     | Sparse marginal forward pass                                                     |
| [`packages/treetime-graph/src/pass.rs`](../../packages/treetime-graph/src/pass.rs)                                                   | Topology-indexed work-first pass storage for Fitch and marginal passes           |
| [`packages/treetime/src/partition/marginal/sparse/message.rs`](../../packages/treetime/src/partition/marginal/sparse/message.rs)     | `combine_messages()` (`#combine_messages`), `propagate_raw()` (`#propagate_raw`) |
| [`packages/treetime/src/seq/mutation.rs`](../../packages/treetime/src/seq/mutation.rs)                                               | Branch mutations by full sequence comparison (`stream_sequence_mutations()`)     |
| [`packages/treetime/src/partition/marginal/sparse/mutations.rs`](../../packages/treetime/src/partition/marginal/sparse/mutations.rs) | Sparse branch mutations (`sparse_edge_mutations()`)                              |
| [`packages/treetime/src/partition/marginal/sparse/partition.rs`](../../packages/treetime/src/partition/marginal/sparse/partition.rs), [`packages/treetime/src/partition/marginal/dense/partition.rs`](../../packages/treetime/src/partition/marginal/dense/partition.rs) | `edge_subs()` of the sparse and dense partitions |
