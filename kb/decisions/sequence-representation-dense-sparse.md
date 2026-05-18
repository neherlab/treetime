# Dense and sparse sequence representation

v1 introduces a sparse sequence representation (`SparseNodePartition` in `packages/treetime/src/partition/sparse.rs:16-84`) alongside the dense representation (`DenseNodePartition` in `packages/treetime/src/partition/dense.rs:16-33`). The sparse path stores probability vectors only for positions where the ancestral state is uncertain, while invariant positions share per-character vectors. This architectural change reduces memory usage proportionally to the fraction of variable sites in the alignment, which ranges from 0.1% to 5% in viral datasets.

v0 uses only dense representation: `seq2array()` (`packages/legacy/treetime/treetime/seq_utils.py:152-204`) converts sequences to NumPy character arrays, and `seq2prof()` (`packages/legacy/treetime/treetime/seq_utils.py:207-229`) converts those to 2D probability matrices of shape `(L, K)` where L is the alignment length and K is the alphabet size (5 for nucleotides including gap, 22 for amino acids). The `SequenceData.make_compressed_alignment()` method (`packages/legacy/treetime/treetime/sequence_data.py:325-464`) groups identical alignment columns and tracks multiplicity, but this column deduplication still stores a full K-width probability vector for each unique column pattern.

The sparse representation affects commands `ancestral`, `timetree`, and `clock`. Users can select dense mode with `--dense=true`; sparse is the default.

## Background: why most alignment positions are invariant

Phylogenetic alignments represent sequences from closely related organisms. Mutations accumulate at a rate proportional to evolutionary time and sequence length. For SARS-CoV-2, the estimated substitution rate is approximately 6.5 x 10^-4 substitutions per site per year [1]. A dataset spanning 5 years of evolution contains roughly 0.3% substitutions per site on average. Most alignment positions remain identical across all sequences.

This sparsity pattern differs from general sparse matrices in scientific computing. Standard sparse formats like CSR (Compressed Sparse Row) store non-zero elements and their indices. Phylogenetic sparsity instead exploits the fact that most positions are invariant (identical character across the tree), allowing those positions to share a single probability vector rather than storing one per position.

Felsenstein's pruning algorithm [2] computes likelihood by traversing the tree postorder (leaves to root), computing partial likelihoods at each node. For a sequence of length L and alphabet size K, each node stores an `(L, K)` matrix of partial likelihoods. The algorithm reduces naive computation from O(4^n) to O(n) for n sequences via dynamic programming. Site pattern compression groups identical columns to avoid redundant computation, but still stores K values per unique pattern at each node.

## v0 implementation

v0 converts sequences to dense NumPy arrays. The `seq2array()` function reads BioPython sequence objects and returns a 1D `np.array` of single characters. For marginal reconstruction, `seq2prof()` converts this to a 2D profile matrix using `profile_map`, a dictionary that maps each character to a one-hot or ambiguity vector. For nucleotides, 'A' maps to `[1,0,0,0,0]`, 'N' maps to `[1,1,1,1,1]`, and IUPAC ambiguity codes map to their corresponding state sets.

Column deduplication in `make_compressed_alignment()` identifies identical alignment columns and tracks their multiplicity. The compressed alignment has L' unique patterns where L' <= L. For a 30,000-position alignment with 1,000 unique column patterns, the compressed length is 1,000. Each unique pattern still carries a full K-width probability vector at every tree node. The multiplicity array weights each pattern's contribution to the total likelihood.

During marginal reconstruction (`treeanc.py:840-932`), the postorder traversal computes `marginal_subtree_LH` at each node: a 2D array of shape `(compressed_L, K)`. Internal nodes multiply propagated child likelihoods. The preorder traversal computes `marginal_outgroup_LH` at each node by combining parent information with sibling subtree likelihoods. All operations use dense matrix arithmetic: `np.dot`, element-wise multiplication, and normalization.

## v1 implementation

v1 stores sequences as `Seq` (`packages/treetime-primitives/src/seq.rs:7-11`), a `Vec<AsciiChar>` with 1 byte per character. This replaces Python's object-per-character overhead with contiguous memory.

### Dense representation

`DenseNodePartition` contains `DenseSeqInfo` (the sequence and gap ranges) and `DenseSeqDistribution` (an `Array2<f64>` of shape `(L, K)`). This matches v0's memory layout. The dense path in `PartitionMarginalDense` (`packages/treetime/src/partition/marginal_dense.rs:24-31`) stores full probability matrices at every node and uses the same belief propagation algorithm as v0.

### Sparse representation

`SparseNodePartition` contains `SparseSeqInfo` (sequence, gap/unknown ranges, Fitch parsimony results) and `SparseSeqDistribution` (`packages/treetime/src/partition/sparse.rs:108-134`). The sparse distribution has:

- `variable: BTreeMap<usize, VarPos>` - position index to probability vector and current state, for positions where parsimony is ambiguous
- `fixed: BTreeMap<AsciiChar, Array1<f64>>` - one shared K-element vector per character type, covering all invariant positions of that character
- `fixed_counts: Composition` - count of invariant positions per character type

The `Composition` struct (`packages/treetime/src/seq/composition.rs:9-12`) tracks how many positions have each character. For an alignment with 25,000 invariant 'A' positions, 3,000 invariant 'C' positions, etc., the fixed contribution to likelihood is computed once per character type and multiplied by its count, rather than computed 25,000 times.

### Compression pipeline

The sparse path runs `compress_sequences()` (`packages/treetime/src/commands/ancestral/fitch.rs:520-546`) before marginal reconstruction. This performs Fitch parsimony (backward and forward passes) to identify variable positions. A position is variable if parsimony cannot unambiguously assign a single state. The result populates `SparseSeqInfo.fitch.variable` with a `BTreeMap<usize, StateSet>` containing only positions that vary across the tree.

Marginal belief propagation (`packages/treetime/src/partition/marginal_passes.rs:16-128`) operates on this compressed representation. The backward pass combines child messages, tracking only variable positions explicitly. Fixed positions contribute through the per-character shared vectors weighted by counts. The forward pass propagates parent information to children, again operating only on variable positions.

### Mutation storage

Sparse edge partitions (`SparseEdgePartition` in `packages/treetime/src/partition/sparse.rs:98-106`) store explicit substitution lists (`Vec<Sub>`) rather than requiring full sequence comparison. This enables efficient extraction of mutations between parent and child nodes.

## Memory comparison

For a tree with N nodes, alignment length L, K alphabet states, and V variable positions:

- v0: `O(N * L' * K)` where L' is the number of unique column patterns
- v1 dense: `O(N * L * K)` (no column deduplication)
- v1 sparse: `O(N * V * K + N * K)` for variable positions plus per-character fixed vectors

With V << L (typically V/L < 0.05 for viral datasets), the sparse representation achieves substantial memory reduction. For a 30,000-position SARS-CoV-2 alignment with ~500 variable sites across 200 nodes and K=5:

- Dense: `200 * 30,000 * 5 * 8 bytes = 240 MB` per partition
- Sparse: `200 * 500 * 5 * 8 + 200 * 5 * 8 = 4 MB + 8 KB` per partition

The sparse representation also reduces computation: belief propagation operations touch V positions rather than L positions.

## Correctness verification

Sparse and dense paths produce identical results. The test suite verifies this in `test_marginal_dense_sparse_log_lh_consistency_gap_free` and `test_marginal_sparse_varpos_matches_dense_profile_gap_free` (`packages/treetime/src/commands/ancestral/__tests__/test_marginal_consistency.rs`). These tests run both paths on the same tree and alignment, asserting that total log-likelihoods match and that variable position probability vectors agree.

## Dense mode selection

The `infer_dense()` function (`packages/treetime/src/partition/algo/infer_dense.rs:1-7`) returns `false` unconditionally. A planned heuristic will select dense mode when tree branches are long enough that most positions become variable, negating the benefit of sparse representation. The threshold depends on branch lengths and mutation rates: when the expected number of mutations per branch approaches the sequence length, sparse representation loses its advantage.

## Practical considerations

- Default behavior uses sparse representation
- Dense mode is available via `--dense=true` for validation or when sparse assumptions break down
- v0's column deduplication is not replicated in v1; the sparse representation provides an alternative memory optimization based on different assumptions
- Sparse representation requires Fitch parsimony as a preprocessing step, adding a tree traversal before marginal reconstruction

## References

[1] SARS-CoV-2 mutation rate estimate from early pandemic sequences. Wikipedia contributors, "Severe acute respiratory syndrome coronavirus 2," Wikipedia. https://en.wikipedia.org/wiki/SARS-CoV-2

[2] Felsenstein J. Evolutionary trees from DNA sequences: a maximum likelihood approach. J Mol Evol. 1981;17(6):368-76. https://doi.org/10.1007/BF01734359. PMID: 7288891.
