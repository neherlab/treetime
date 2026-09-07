# Mutation-first sequence representation

## Background

Every backend materializes a full sequence for every node. `SparseSeqInfo::sequence` ([packages/treetime/src/partition/storage/sparse.rs](../../packages/treetime/src/partition/storage/sparse.rs)) holds one `Seq` per node, and after [the parsimony-chain change](../decisions/ancestral-marginal-sparse-parsimony-chain.md) tips additionally hold their emitted sequence in `SparseNodePartition::emitted`. Both scale as `O(N * L)` in node count and alignment length.

That is redundant. A node's sequence is fully determined by the root sequence plus the per-edge substitutions and indels on the path to it, plus the node's own missing-data mask. All of those are already stored, and all of them are compact: substitutions live in `SparseEdgePartition::subs_ml`/`subs_fitch`, indels in `SparseEdgePartition::indels`, and missing data in the range tracks `SparseSeqInfo::{unknown, gaps, non_char}`.

The sparse backend exists to avoid `O(N * L)` storage of probability vectors. It then spends `O(N * L)` on sequences anyway.

## Evidence

Measured on `data/sc2/4500` (10199 nodes, L = 29903, mean 2.6 substitutions per edge):

| Component | Size |
| --------- | ---- |
| Per-node sequences | 305 MB |
| `emitted` for tips | 153 MB |
| Root sequence | 0.03 MB |
| All edge substitutions (26300, 8 B each) | 0.21 MB |
| All unknown + gap ranges (121980, 16 B each) | 1.95 MB |
| Tip IUPAC positions (2510, 9 B each) | 0.02 MB |

Roughly 458 MB against 2.2 MB, a factor of ~200. The ratio grows with alignment length at fixed divergence, which is the regime the sparse backend targets.

### The reported mutations are sufficient

This was checked rather than assumed. For each of the 10198 edges, the emitted parent and child sequences were compared position by position and the differences matched against the substitutions reported on that edge:

- **0** edges carry a reported substitution that is not an actual difference.
- **0** positions differ between parent and child, with both endpoints canonical, without a reported substitution.

Every unreported difference has `N`, `-`, or an IUPAC code at one end: maskings (`T -> N`, `A -> -`), insertions (`- -> T`), and observed tip ambiguity (`G -> R`). Each of those is carried by a track that is already stored separately.

This follows from how the ranges are built and is not an accident of the dataset. An internal node's `non_char` is the intersection of its children's, so a parent cannot be masked at a position where a child holds a residue; the "parent masked, child canonical" case that would lose a substitution is unreachable.

A naive check that walks `root + path substitutions` and compares against emitted sequences appears to show thousands of canonical mismatches at internal nodes. That is accumulated drift, not information loss: once a masked position is left unapplied, every descendant inherits the stale value and compares canonical against canonical downstream. The per-edge comparison is the correct diagnostic.

### Reconstruction rule

```
sequence(node) = root_sequence
               + substitutions along the root -> node path
               + indels along the path
               + node's unknown ranges          (N mask)
               + node's gap ranges              (deletion mask)
               + observed IUPAC codes           (tips only, seq.fitch.variable)
```

## Proposed change

Stop storing per-node sequences. Keep the root sequence, the per-edge mutation and indel tracks, and the per-node range tracks. Materialize a sequence only where one is genuinely required, and do it transiently.

The natural mechanism is a preorder traversal carrying a single mutable buffer: apply the edge's substitutions and indels on the way down, undo them on the way back up. Working memory becomes `O(L)` rather than `O(N * L)`. `ancestral_reconstruction_marginal` already traverses with `iter_depth_first_preorder_forward`, which is a genuine DFS and single-threaded.

## Design axes

### Axis 1: the random-access accessor

`PartitionBranchOps::node_sequence(&self, node_key) -> Seq` is called from about ten sites (`commands/shared/tree_output.rs`, `commands/ancestral/augur_node_data.rs`, `partition/io/augur.rs`, `partition/timetree/branch.rs`) in arbitrary key order. It is the central obstacle: a transient buffer cannot serve arbitrary-order queries.

- **Option A - streaming visitor.** Replace the accessor with a traversal that hands each node its sequence in preorder. Matches how the output layer already walks the tree. Largest diff, cleanest end state.
- **Option B - on-demand replay.** Keep the accessor and reconstruct by walking root-to-node applying diffs, `O(depth * muts)` per call. Small diff, but per-node queries across the whole tree make it quadratic-ish.
- **Option C - per-node diff map.** Store `BTreeMap<usize, AsciiChar>` of the diff from the root per node. Preserves random access at `O(depth * 2.6)` entries per node instead of `L`. Still a large win, much smaller change than A.

### Axis 2: sampling

`--sample-from-profile=all` redraws every position from its fixed block ([reconstruct.rs](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs)), so a realization is a full-length random object that is not expressible as a small diff. Options: restrict that mode to streaming output only; or keep a materialized sequence for sampled nodes, paying `O(N * L)` only in that mode. `argmax` (the default) and `root` are unaffected.

### Axis 3: the augur `sequence` field

`AugurNodeDataJsonAncestralNode::sequence` is already `Option<String>` in the schema. Augur pipelines consume `muts` plus the root `reference`. Gating the per-node sequence behind a flag removes the largest consumer without changing the default output contract. Needs a compatibility check against `augur export v2`.

### Axis 4: tips

Tips are half the nodes here (5100 of 10199) and their sequences are the input alignment held a second time. They could reference the alignment rather than copy it, or be diff-encoded like internal nodes.

### Axis 5: two mutation sets

If a sequence is *defined* as root plus path, `parent + mutations == child` becomes a tautology and [`test_marginal_tip_parent_plus_muts_equals_child`](../../packages/treetime/src/ancestral/__tests__/test_marginal_tip_reconstruction.rs) loses its content. The reported set (`subs_ml`, filtered for augur parity) and the reconstruction set are then distinct concepts and should be distinct in the types, or a future reader will reconstruct from the wrong one.

## Expected impact

- Peak memory on `data/sc2/4500` drops by roughly 450 MB for the ancestral path.
- `capture_ancestral_states` ([timetree/convergence/sequence_changes.rs](../../packages/treetime/src/timetree/convergence/sequence_changes.rs)) snapshots every internal sequence twice per timetree iteration, about 610 MB on this dataset, purely to count changed positions. Comparing per-node mutation sets is both cheaper and a more direct convergence signal.
- Dense is unaffected; it stores `(L, K)` profiles per node regardless, and per-node sequences are not its dominant term.
- No output changes if the streaming traversal produces byte-identical sequences, which is the acceptance criterion below.

## Suggested staging

Each step stands alone and is independently revertable:

1. **Convergence check** - compare mutation sets rather than sequences. Self-contained, no output change, immediate win.
2. **Gate the augur per-node `sequence`** behind a flag, after checking `augur export v2` compatibility.
3. **Streaming emit traversal** for the reconstructed FASTA and the tree-annotation outputs.
4. **Drop `seq.sequence` and `emitted`** once nothing random-accesses them, resolving axes 1, 2 and 4.

## Validation plan

- Byte-identical reconstructed FASTA, augur node-data JSON, auspice, nexus and PhyloXML against the current implementation on `wrong_reconstruction/subset` (676 tips, 29903 columns) and `data/sc2/4500`.
- Byte-identical output under `--impute-missing-data`, `--sample-from-profile root --seed 42` and `--sample-from-profile all --seed 42`.
- Sparse still byte-identical to dense on the four `wrong_reconstruction` alignments.
- `optimize` and `timetree` tree and node-data unchanged, except where step 1 deliberately changes the convergence metric, which should be reported as an iteration-count difference rather than an output difference.
- Peak RSS measured before and after on `data/sc2/4500`.

## Related issues

- [M-ancestral-sparse-stores-a-full-sequence-per-node](../issues/M-ancestral-sparse-stores-a-full-sequence-per-node.md)

## Open questions

- Does any downstream consumer (augur, nextstrain) rely on the per-node `sequence` field in node-data JSON?
- Should the parsimony chain and the MAP sequence remain distinguishable at the API level, or is only the MAP sequence ever requested outside the marginal passes?
- Is `--sample-from-profile=all` used enough to justify keeping a materialization path for it?
