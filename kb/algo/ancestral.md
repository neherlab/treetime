# Ancestral Reconstruction Algorithms

## Fitch Parsimony

Maximum parsimony (<a id="cite-1"></a>[Fitch 1971](https://doi.org/10.2307/2412116) [[1](#ref-1)]) reconstructs ancestral character states by minimizing the total number of state changes on the tree. The method makes no assumptions about branch lengths or substitution rates, treating all state transitions as equally costly. It remains widely used for seeding ML optimization with initial ancestral assignments due to its speed, and for compression of sequence data to variable-position-only representations.

v1: `fn fitch_backward()`, `fn fitch_forward()`, and `fn fitch_cleanup()` in [`packages/treetime/src/partition/fitch/passes.rs#L90-L305`](../../packages/treetime/src/partition/fitch/passes.rs#L90-L305), driven per partition by [`packages/treetime/src/ancestral/fitch.rs`](../../packages/treetime/src/ancestral/fitch.rs).
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

- `fitch_backward()` (`#fitch_backward`) and `FitchBackwardPass::node()` (`#FitchBackwardPass::node`): postorder pass computing state sets and parsimony score
- `fitch_forward()` (`#fitch_forward`) and `fitch_forward_node()` (`#fitch_forward_node`): preorder pass resolving ambiguities

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

Both passes run on the graph scans of [`packages/treetime-graph/src/schedule_scan.rs`](../../packages/treetime-graph/src/schedule_scan.rs): the backward pass on `GraphSchedule::scan_up`, the forward pass on `GraphSchedule::scan_down`, after `require_tree()` checks that the graph is a tree. Each scan takes the per-node inputs by value, gives each node step only its own input and the messages of its children (up) or its parent (down), and returns the per-node outputs. A node runs as soon as its dependencies complete rather than at a per-level barrier, and children are combined in edge order, so results do not depend on the thread count.

### v1 implementations

**Dense** (all positions): `fn indexed_backward()` and `fn indexed_forward()` in [`packages/treetime/src/partition/marginal/shared/pass.rs`](../../packages/treetime/src/partition/marginal/shared/pass.rs) (shared with discrete partitions), partition in [`packages/treetime/src/partition/marginal/dense/partition.rs`](../../packages/treetime/src/partition/marginal/dense/partition.rs). Stores full probability vectors at every alignment position. Used when the full profile is needed (e.g., GTR inference from data).

**Sparse** (variable positions only): [`packages/treetime/src/partition/marginal/sparse/`](../../packages/treetime/src/partition/marginal/sparse/) (`backward.rs`, `forward.rs`, `message.rs`). Stores explicit profiles only at positions that vary near each node; see [Sparse marginal design](#sparse-marginal-design). Much faster for conserved alignments where >90% of positions are invariant.

v0: [`packages/legacy/treetime/treetime/treeanc.py#L762-L927`](../../packages/legacy/treetime/treetime/treeanc.py#L762-L927).

### Sparse node data model

Each sparse node has observations (`struct SparseNodeObs`) and a pass state (`struct SparseNodeState`), both in [`packages/treetime/src/partition/storage/sparse.rs`](../../packages/treetime/src/partition/storage/sparse.rs):

Observations: `composition` (character counts), `gaps`/`unknown`/`non_char` (position ranges where the node has no informative data), and the Fitch state. MAP reconstruction replaces the node sequence and recomputes `composition` atomically, so later sparse passes use counts from the current reconstruction rather than the earlier Fitch assignment.

Fitch state: `fitch.variable` (ambiguous positions from parsimony, as state sets), `fitch.chosen_state` (resolved state per variable position from Fitch forward pass). After the Fitch passes, `fn fitch_cleanup()` [`packages/treetime/src/partition/fitch/passes.rs#L298-L305`](../../packages/treetime/src/partition/fitch/passes.rs#L298-L305) clears `fitch.variable` for internal nodes, so only leaves keep it.

Pass state: `sequence` (stored sequence) and `profile` (`struct SparseSeqDistribution`). A profile has two parts. `variable` (`struct VarPosMap`, [`packages/treetime/src/partition/storage/var_pos_map.rs`](../../packages/treetime/src/partition/storage/var_pos_map.rs)) lists the explicit positions sorted by position; each entry holds the Fitch `state` of the position and its probability vector `dis`, stored inline in the entry, and a lookup returns a `VarPos { dis, state }` view of it. `fixed` (`struct FixedProfiles`) lists the fixed states sorted by state; each entry holds one shared vector, used for every position of that state not listed in `variable`, and the `count` of those positions. A position that is certainly `A` in every message needs no vector of its own; it adds one to the count of `A` and uses the vector of `A`. Vectors of up to 4 states (`INLINE_STATES`, the nucleotide alphabet) need no heap allocation; see [kb/issues/N-ancestral-sparse-profiles-allocate-for-large-alphabets.md](../issues/N-ancestral-sparse-profiles-allocate-for-large-alphabets.md).

`composition` gives the fixed counts of `msg_to_child` (`fn compute_msg_to_child()`): the count of each parent state minus the explicit positions with that parent state and minus the child-missing sites with that parent state (see [Sparse missing data](#sparse-missing-data)). These counts are the per-character multiplicity weights for branch-length optimization ([`packages/treetime/src/partition/optimize/sparse.rs`](../../packages/treetime/src/partition/optimize/sparse.rs)). `non_char` excludes positions from `edge_effective_length` and from the marginal substitution estimates of an edge.

The `non_char` of an internal node mixes two kinds of positions: positions where its subtree has no data, and gaps that the Fitch indel pass resolves although a child subtree has data there (`fn resolve_indels_backward()` in [`packages/treetime/src/seq/indel.rs`](../../packages/treetime/src/seq/indel.rs) resolves a gap when every child is gapped, unknown, or variable in indels). The marginal passes therefore do not read `non_char` to decide whether a subtree has data; they read the no-data range of [Sparse missing data](#sparse-missing-data).

Non-char (N, gap) differences between parent and child are tracked through `non_char` ranges, not as Fitch substitutions. Applying a child edge's fitch_subs and indels to the parent's composition does not reproduce the child's composition when non-char positions differ. This is by design: Fitch parsimony assigns concrete states at ambiguous positions, and N/gap masking is a separate layer.

### Sparse marginal design

Branches in typical TreeTime datasets are short, and on a given edge or near a given node only a small minority of sites change. Computing propagators for every site wastes most of the work on sites with negligible probability of change. Two alternatives exist. Parsimony avoids probabilistic models altogether. MAPLE (<a id="cite-8"></a>[De Maio et al. 2023](https://doi.org/10.1038/s41588-023-01368-0) [[8](#ref-8)]) stores explicit vectors only for positions that differ from a reference genome, but that set grows with the distance from the reference, so the representation stops being sparse as the tree gets deeper. v1 instead keeps a site explicit only where it varies near the current node: at Fitch substitutions on adjacent edges, at ambiguous leaf sites, and where the messages disagree or are uncertain. Every other site is represented by the shared `fixed` vector of its Fitch state. The Fitch reconstruction supplies everything the passes need: the root sequence, the substitutions on each edge, the gap and unknown ranges of each node, the ambiguous sites of each leaf with their chosen states, and the state counts of each node.

**Node sequences** (`fn fitch_node_sequences()` [`packages/treetime/src/partition/marginal/sparse/sequences.rs#L12-L50`](../../packages/treetime/src/partition/marginal/sparse/sequences.rs#L12-L50)): every sparse marginal update starts with a preorder pass that rebuilds the Fitch sequence of each internal node from its parent's sequence with the Fitch substitutions and indels of its edge and its own unknown ranges (`fn parsimony_seq()`). The root takes `root_sequence`, and a leaf keeps its stored sequence. `fn into_marginal_sparse()` stores sequences only for the root and the leaves, and the nodes that a merge or a hoist of `optimize` creates start empty, so this pass gives every node its Fitch state before the backward pass reads it. The pass reuses a stored sequence that equals the rebuilt one, so memory stays at one sequence per node.

**Backward pass** (`fn SparseBackwardPass::node()` [`packages/treetime/src/partition/marginal/sparse/backward.rs#L66-L116`](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L66-L116)):

- Leaf (`fn leaf_message()`): `fixed` holds the profile of each determined state (canonical and ambiguous) with the leaf composition count of that state. Each ambiguous position becomes an explicit position whose `dis` is the profile of its ambiguity set and whose `state` is the Fitch chosen state (or one member of the set when the chosen state is not canonical)
- Internal node, collection (`fn combine_child_messages()` [`packages/treetime/src/partition/marginal/sparse/backward.rs#L133-L206`](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L133-L206)): the explicit positions of the node are (1) every gapped-parent site of the node ([Sparse missing data](#sparse-missing-data)), (2) every position with a Fitch substitution on a child edge, with the substitution's reference state as the node state and its query state recorded for that child, (3) every position explicit in a child message, with the child's state, and (4) every child-missing site of a child edge where the node's Fitch state is canonical, with that state. Substitutions take precedence. The node state of a position in the node's `non_char` is then the gap (`fn own_state()`)
- Internal node, child states: for each child and each collected position without a substitution on that child's edge, the child state is the child's own Fitch state where it is canonical and the node state otherwise, except in the child's no-data range, where it is gap or unknown. The child state selects the `fixed` vector used when the position is not explicit in that child's message. This is why `combine_messages()` takes a reference state per message rather than per child: the same function serves the forward pass, where the two messages come from the parent and the child
- Candidate site table: both passes and `compute_msg_to_child()` record the collection as events in one list (`struct CandidateSiteEvents` in [`packages/treetime/src/partition/marginal/sparse/candidate_sites.rs`](../../packages/treetime/src/partition/marginal/sparse/candidate_sites.rs)). Each event sets the node state or the reference state of one message at a position, either replacing an earlier value or only filling an empty slot, as a map insert or `entry().or_insert()` would. A stable sort by position keeps the recording order within a position, and `finish()` resolves each slot and emits one row per position that has a node state (`struct CandidateSites`)
- `combine_messages()` then multiplies the messages and propagates the result to the parent edge with $P(t)^T$. At the root it also multiplies in the equilibrium frequencies $\pi$
- Missing data: a gap or unknown reference adds no factor, so a child contributes nothing at its child-missing sites and at gapped-parent sites where it has no data. Such a site has a gap or unknown reference, so it stays explicit at the node and folds into the `fixed` vector of its state at an ancestor only under the rule below. A gapped-parent site has no canonical node state: it stays explicit, labeled with its most likely state, and leaves the per-state weights unchanged, because `composition` does not count it. The root adds these sites to the tree likelihood in the same way

**Message combination** (`fn combine_messages()` [`packages/treetime/src/partition/marginal/sparse/message.rs#L22-L107`](../../packages/treetime/src/partition/marginal/sparse/message.rs#L22-L107)): it takes the candidate site table, with one reference state per incoming message. For each site, it sums the log of the message vectors in log space: the explicit vector of a message where the position is explicit in it, otherwise the `fixed` vector of the message's reference state (a gap or unknown reference contributes nothing), otherwise the `fixed` vector of the site state. It normalizes the sum in place with `softmax_with_log_norm_inplace()`. A position stays explicit unless both conditions hold:

- the largest probability exceeds $1 - \epsilon$ with $\epsilon = 10^{-4}$ (`fn is_site_resolved()`)
- every message agrees with the node's Fitch state at the position

The check does not require the argmax to equal the Fitch state; both conditions together make a disagreement unlikely but not impossible. A dropped position counts toward the `fixed` count of its state. The `fixed` vectors are combined once per canonical state. The log-likelihood weights each by the composition count of the state minus every candidate site of that state, while the stored `count` subtracts only the sites that stay explicit. Dropping positions is an approximation; see [kb/decisions/ancestral-marginal-sparse-parsimony-chain.md](../decisions/ancestral-marginal-sparse-parsimony-chain.md) and [kb/issues/H-ancestral-sparse-dense-internal-residues-diverge.md](../issues/H-ancestral-sparse-dense-internal-residues-diverge.md).

**Forward pass** (`fn SparseForwardPass::node()` and `fn SparseForwardPass::edge()` [`packages/treetime/src/partition/marginal/sparse/forward.rs#L86-L202`](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L86-L202)):

- `msg_to_child` (`fn compute_msg_to_child()` [`packages/treetime/src/partition/marginal/sparse/forward.rs#L221-L310`](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L221-L310)) divides the parent profile by the child's `msg_from_child`, position by position and for each `fixed` vector, and normalizes. Its explicit positions are the edge substitutions and the explicit positions of the parent profile and of `msg_from_child`, without the child's no-data range
- Child-missing sites in `msg_to_child`: every position of the child's no-data range that is not explicit is removed from the `fixed` count of the parent's state, without a division and without a normalization term, because the child contributed nothing there and the outside message is the parent posterior with normalization 1. The branch-length objective reads these counts, so the sites leave it
- `msg_to_child` is propagated along the edge with $P(t)$ to give `msg_from_parent`
- The node profile is `combine_messages(msg_from_parent, msg_to_parent)` without a root prior. Its explicit positions are the edge substitutions (child state), the explicit positions of `msg_to_child` outside the child's no-data range, and the explicit positions of `msg_to_parent`
- Outside messages at no-data sites (`PartitionMarginalSparse::impute_missing_data`): when a run imputes missing data, the pass also keeps the undivided parent posterior at every no-data site of the child, propagated across the edge: in `msg_from_parent` for a leaf, where `fn leaf_posterior()` reads it to impute the leaf's `N` sites, and in the node profile for an internal node, so the message passes down a clade without data to its leaves. `ancestral` and `timetree` set the flag when they impute, `optimize` never sets it, because these per-site entries serve only imputation
- `fn compute_ml_subs_for_nodes()` places no marginal substitution at a position in the `non_char` of either node or in the child's no-data range

Both passes also handle per-site rates (`fn propagate_raw_per_site()`) and the log-likelihood bookkeeping of the division step.

### Sparse missing data

The **no-data range** of a node is the set of positions where its subtree has no data: the gaps and unknown positions of a leaf, and the intersection of the children's no-data ranges at an internal node. The backward pass computes it bottom-up in every update and stores it on the edge to the parent (`SparseEdgeBackward.child_no_data`), so it is never stale after a topology change. Two kinds of sites follow from it:

- **Child-missing site** of an edge: a position in the child's no-data range and outside the parent's `non_char`. The child's subtree has no data there, while the parent has data from other children. The child must contribute nothing at such a site in the backward pass, the forward pass, and the branch-length objective; the shared `fixed` vector of the parent's state would otherwise count it as an unchanged site and pull the branch length toward zero
- **Gapped-parent site** of a node: a position in the node's `non_char` and outside its no-data range, a Fitch-resolved gap with data below. The data of the children with data there must enter the likelihood, although `composition` does not count the site

The dense backend treats gaps and `N` as no data at every node, as v0 does, and is the oracle for both cases: `packages/treetime/src/partition/marginal/sparse/__tests__/test_missing_data.rs` compares the root log-likelihood, the change of every edge objective and of the tree log-likelihood with the branch length, the node log-likelihoods, and the leaf posteriors at `N` sites with dense on cases with `N` and gap runs in leaves and clades, an all-`N` leaf, and gapped ancestors below the root and at the root.

### v0 differences

v1 backward pass uses log-space arithmetic with logsumexp normalization (dense `normalize_from_log()`, sparse `softmax_with_log_norm_inplace()` in `combine_messages()`). v1 forward pass uses plain probability space (division). v0 uses neg-log space throughout. v1 dense uses deterministic `argmax_first()` (`#argmax_first`) for sequence extraction (leftmost state wins ties); v0 uses `np.argmax()` which has undefined tie-breaking.

### Key functions

- `SparseBackwardPass::node()` (`#SparseBackwardPass::node`, sparse) and `IndexedBackwardPass::node()` (`#IndexedBackwardPass::node`, dense): compute partial likelihoods via GTR matrix multiplication
- `SparseForwardPass::edge()` (`#SparseForwardPass::edge`, sparse) and `IndexedForwardPass::node()` (`#IndexedForwardPass::node`, dense): compute outgroup messages via cavity/division
- `combine_messages()` (`#combine_messages`): combines child messages in log-space via logsumexp normalization
- `propagate_raw()` (`#propagate_raw`): GTR matrix-vector product `P(t) * profile`

### Complexity

O(n _ k^2 _ L) total for n nodes, k alphabet states (4 for nucleotides, 20 for amino acids), and L alignment positions. Brute force summation over all possible ancestral assignments would be O(k^n \* L), exponential in tree size.

v0 reduces L to the number of distinct alignment columns by grouping identical columns. The grouping pays only when columns repeat; on alignments where almost every column differs, its bookkeeping costs more than it saves. v1 dense does not group columns, and v1 sparse avoids per-column work for invariant positions instead: [kb/decisions/sequence-representation-dense-sparse.md](../decisions/sequence-representation-dense-sparse.md).

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

| File                                                                                                                                                                                                                                                                     | Algorithms                                                                       |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------- |
| [`packages/treetime/src/partition/fitch/passes.rs`](../../packages/treetime/src/partition/fitch/passes.rs)                                                                                                                                                               | Fitch parsimony (backward, forward, cleanup)                                     |
| [`packages/treetime/src/ancestral/fitch.rs`](../../packages/treetime/src/ancestral/fitch.rs)                                                                                                                                                                             | Fitch reconstruction over all partitions                                         |
| [`packages/treetime/src/ancestral/pipeline.rs`](../../packages/treetime/src/ancestral/pipeline.rs)                                                                                                                                                                       | Ancestral pipeline, method dispatch                                              |
| [`packages/app-commands/src/commands/ancestral/run.rs`](../../packages/app-commands/src/commands/ancestral/run.rs)                                                                                                                                                       | Ancestral command entry point                                                    |
| [`packages/treetime/src/partition/marginal/shared/pass.rs`](../../packages/treetime/src/partition/marginal/shared/pass.rs)                                                                                                                                               | Dense and discrete marginal passes (Felsenstein pruning)                         |
| [`packages/treetime/src/partition/marginal/sparse/backward.rs`](../../packages/treetime/src/partition/marginal/sparse/backward.rs)                                                                                                                                       | Sparse marginal backward pass                                                    |
| [`packages/treetime/src/partition/marginal/sparse/forward.rs`](../../packages/treetime/src/partition/marginal/sparse/forward.rs)                                                                                                                                         | Sparse marginal forward pass                                                     |
| [`packages/treetime/src/partition/marginal/sparse/sequences.rs`](../../packages/treetime/src/partition/marginal/sparse/sequences.rs)                                                                                                                                     | Node sequences before the sparse backward pass                                   |
| [`packages/treetime-graph/src/schedule_scan.rs`](../../packages/treetime-graph/src/schedule_scan.rs)                                                                                                                                                                     | Parallel graph scans that run the Fitch and marginal passes                      |
| [`packages/treetime/src/partition/marginal/sparse/message.rs`](../../packages/treetime/src/partition/marginal/sparse/message.rs)                                                                                                                                         | `combine_messages()` (`#combine_messages`), `propagate_raw()` (`#propagate_raw`) |
| [`packages/treetime/src/seq/mutation.rs`](../../packages/treetime/src/seq/mutation.rs)                                                                                                                                                                                   | Branch mutations by full sequence comparison (`stream_sequence_mutations()`)     |
| [`packages/treetime/src/partition/marginal/sparse/mutations.rs`](../../packages/treetime/src/partition/marginal/sparse/mutations.rs)                                                                                                                                     | Sparse branch mutations (`sparse_edge_mutations()`)                              |
| [`packages/treetime/src/partition/marginal/sparse/partition.rs`](../../packages/treetime/src/partition/marginal/sparse/partition.rs), [`packages/treetime/src/partition/marginal/dense/partition.rs`](../../packages/treetime/src/partition/marginal/dense/partition.rs) | `edge_subs()` of the sparse and dense partitions                                 |
