# Dense and sparse sequence representation

| Property    | Value                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Intentional architectural change                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| v1 Location | `DenseNodePartition` (`#DenseNodePartition`) in [`packages/treetime/src/representation/payload/dense.rs`](../../packages/treetime/src/representation/payload/dense.rs), `SparseNodePartition` (`#SparseNodePartition`) in [`packages/treetime/src/representation/payload/sparse.rs`](../../packages/treetime/src/representation/payload/sparse.rs), `Seq` (`#Seq`) in [`packages/treetime-primitives/src/seq.rs`](../../packages/treetime-primitives/src/seq.rs) |
| v0 Location | `seq2array()` (`#seq2array`) in [`packages/legacy/treetime/treetime/seq_utils.py`](../../packages/legacy/treetime/treetime/seq_utils.py), `SequenceData` (`#SequenceData`) in [`packages/legacy/treetime/treetime/sequence_data.py`](../../packages/legacy/treetime/treetime/sequence_data.py)                                                                                                                                                                   |
| Affects     | Memory usage, runtime performance on large alignments                                                                                                                                                                                                                                                                                                                                                                                                            |
| Commands    | `ancestral`, `timetree`, `clock`                                                                                                                                                                                                                                                                                                                                                                                                                                 |

## What v0 does

v0 stores all sequence data as dense NumPy arrays. Leaf sequences are converted to `np.array` of single characters via `seq2array()`. During marginal reconstruction, each node carries a `marginal_subtree_LH` profile: a 2D `np.array` of shape `(L, K)` where `L` is the (compressed) alignment length and `K` is the alphabet size (5 for nucleotides, 22 for amino acids). Every position in the alignment has a full probability vector, regardless of whether that position is variable across the tree.

v0 reduces computational cost through column deduplication (`SequenceData.make_compressed_alignment()`): identical alignment columns are merged and tracked with a multiplicity count. This reduces `L` when many columns share the same pattern, but each remaining column still carries a full `K`-wide probability vector at every node.

All belief propagation operations (postorder, preorder) use dense matrix arithmetic: `np.dot`, element-wise multiplication, and normalization over full `(L, K)` matrices.

## What v1 does

v1 provides two representations, selectable via the `--dense` CLI flag (default: sparse).

**Dense** (`PartitionMarginalDense`): stores a full `Array2<f64>` of shape `(L, K)` per node in `DenseSeqDis.dis`, matching v0's memory layout. Belief propagation operates on full matrices. Sequences are stored as `Seq` (`Vec<AsciiChar>`) in `DenseSeqInfo`.

**Sparse** (`PartitionMarginalSparse`): stores probability vectors only for variable positions. Each node carries a `MarginalSparseSeqDistribution` with:

- `variable: BTreeMap<usize, VarPos>` - per-position probability vector and current state for positions where the ancestral state is uncertain
- `fixed: BTreeMap<AsciiChar, Array1<f64>>` - one shared probability vector per character type, covering all invariant positions of that type
- `fixed_counts: Composition` - count of invariant positions per character

The sparse path first runs Fitch parsimony to identify variable positions and compute a compressed representation (`compress_sequences()`), then performs marginal belief propagation only on variable positions plus per-character summaries for fixed positions.

Sequences at all nodes are stored as `Seq` (`Vec<AsciiChar>`) in `SparseSeqInfo`, with mutations represented as `Vec<Sub>` on edges rather than full sequences at every node.

## Why v1 changes this

Phylogenetic alignments are dominated by invariant positions. In a typical viral dataset, fewer than 5% of alignment columns are variable. v0's dense approach allocates a full `(L, K)` matrix at every node regardless, leading to memory usage proportional to `L * K * N_nodes`.

The sparse representation reduces memory to roughly proportional to `V * K * N_nodes + L`, where `V` is the number of variable positions (V << L). For a 30,000-position SARS-CoV-2 alignment with ~500 variable sites across 200 nodes, this is approximately a 60x reduction in profile storage.

The per-character fixed-position summary (`fixed` and `fixed_counts`) preserves the log-likelihood contribution of invariant sites without storing per-position vectors. Belief propagation at fixed positions reduces to a single matrix-vector product per character type rather than per position.

v1 retains the dense path for correctness validation and for cases where sparsity assumptions break down (long branches producing many variable positions).

## Practical impact

- Default behavior uses sparse representation. Users can force dense mode with `--dense=true`.
- Sparse and dense produce identical log-likelihoods and ancestral sequences (verified by `test_marginal_dense_sparse_log_lh_consistency_gap_free` and `test_marginal_sparse_varpos_matches_dense_profile_gap_free` in [`packages/treetime/src/commands/ancestral/__tests__/test_marginal_consistency.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_consistency.rs)).
- v0's column deduplication compression is replaced by Fitch-based variable position identification in the sparse path.
- Edge data in sparse mode stores explicit substitution lists (`Vec<Sub>`) rather than requiring full sequence comparison to extract mutations.
- The `infer_dense()` function in [`packages/treetime/src/representation/algo/infer_dense.rs`](../../packages/treetime/src/representation/algo/infer_dense.rs) currently returns `false` unconditionally. A planned heuristic will select dense mode when tree branches are long enough that most positions become variable.
