# The sparse stored sequence is the parsimony chain

This document records two correctness fixes in the v1 sparse marginal ancestral reconstruction and the representation choice that resolves both. They affect internal-node sequences emitted by `treetime ancestral --method-anc=marginal` and `treetime timetree` under the default sparse backend. Posteriors, log-likelihoods, and the dense backend are unchanged.

This is the internal-node counterpart of the tip corruption recorded in [Marginal Tip Reconstruction and Missing-Data Imputation](ancestral-marginal-tip-reconstruction-and-imputation.md). Both come from the same source - reconstructing a node by chaining off its parent - but the tip fix removed leaves from the chain and left the internal-node case in place.

## Background

The sparse backend does not store a posterior at every position. `combine_messages` (`packages/treetime/src/partition/marginal/sparse/message.rs:17`) keeps an explicit distribution only for positions that stay genuinely uncertain; once a position's posterior concentrates on one state above `1 - EPS` (`EPS = 1e-4`, `message.rs:15`) and every incoming message agrees with the parsimony reference state, the position is dropped from `profile.variable` and represented by the shared per-character `fixed` vector instead.

Reconstruction used to rebuild each node's sequence by cloning the parent's _reconstructed_ sequence, applying the parsimony substitutions and indels on the edge, and overwriting the positions the node still held in `profile.variable`. Positions absent from `variable` were taken from the parent unchanged.

## Problem: a deviating state propagates through a whole clade

The chain assumed the parent's own state equals its parsimony reference state. Marginal reconstruction is free to disagree with parsimony and does so systematically around near-zero branch lengths: a tip on a zero-length branch is nearly conclusive evidence about its parent, so the maximum-likelihood explanation is that the parent already carried the tip's state, while parsimony places the substitution on the tip branch.

At such a node the reconstruction correctly wrote the deviating state, because the position was still in its `profile.variable`. But every descendant that had resolved the position inherited that state instead of its own. One node's legitimate state change became a subtree-wide one, and the deeper and more conserved the clade, the more nodes it flipped - resolving a position out of `variable` is exactly what a deep, conserved clade does.

### Worked example

On a 676-tip SARS-CoV-2 subtree with a 1000-column alignment starting at genome position 22000, position 928 holds `C` in 674 tips and `T` in two (`JN.1.11`, `JN.1.11.2`). Tip `JN.1.11` hangs off `NODE_0000080` on a zero-length branch. Fitch places the single `C928T` on that tip branch, but the near-zero branch makes `T` at the parent more likely than a mutation: the true posterior at `NODE_0000080` is `C = 0.4680 / T = 0.5320`.

The sparse backend computed that posterior exactly - verified against an independent Felsenstein implementation and against the dense backend, all three agreeing to every printed digit. `NODE_0000080` and `NODE_0000081` are genuinely `T`. The clade below them is all-`C` with upward evidence of order `1e-32`, so its nodes resolved position 928 out of `variable` and inherited `T` from above. The run emitted 194 internal `T` where only 2 are correct.

### Why it looked length-dependent

The apparent dependence on alignment length is indirect. The default `--model infer` estimates the GTR from the data, so adding alignment columns shifts the rate matrix, which moves a near-tied posterior across the 50/50 line. On this dataset the flip is caused by a single column: truncating to 985 columns reconstructs position 928 correctly, and adding column 986 (itself a `C`/`T` polymorphism) flips `NODE_0000080` and with it the whole clade. `--model jc69` never flips it. Longer alignments make the flip more likely and can move it higher in the tree, where the corrupted clade is larger.

Zero-length branches are clamped to `MIN_BRANCH_LENGTH_FRACTION / seq_length` (`packages/treetime/src/constants.rs:7`, `packages/treetime/src/hacks/fix_branch_length.rs`), so the evidence a zero-branch tip contributes about its parent grows with alignment length. That affects both backends identically and is not part of this fix, but it is what makes these near-tied parent states common enough to matter.

## Problem: inherited deletions were overwritten

The same routine re-applied deletions _last_, to stop a posterior from resurrecting a residue at a deleted site - but only the deletions recorded on the node's own parent edge. A node whose gap was inherited from further up carries no indel on its own edge, so nothing masked the position again and the posterior's residue stood.

On `data/sc2/4500` position 28369 (4332 of 5100 tips deleted), the old reconstruction emitted `A` at 5014 internal nodes. Dense reports 4555 gaps and 544 `A`. At position 23009 it emitted `A` at 864 nodes that dense reports as gaps. Both are now fixed.

Parsimony reconstruction and `--dense=true` were never affected by either defect: neither chains a node's sequence off its parent this way.

## Decision

`SparseSeqInfo::sequence` always holds the **parsimony sequence**: the root sequence with each edge's Fitch substitutions and indels applied, masked by the node's own missing data. `parsimony_seq` (`packages/treetime/src/partition/marginal/sparse/reconstruct.rs:23`) builds it once in the marginal forward pass, and nothing rewrites it afterwards.

The MAP sequence is a derived view: `map_seq` (`reconstruct.rs:50`) takes the parsimony sequence and resolves each position in `profile.variable` from its posterior. Positions absent from `variable` resolved to a single state during message combination, and that state is by construction the parsimony state, so taking them from the chain is exact.

Both defects follow from the separation rather than being patched:

- The chain never carries a MAP state, so a node whose argmax disagrees with Fitch has nothing to propagate.
- A variable position whose parsimony character is already a gap is skipped, which covers the node's full inherited gap set rather than just its own edge's deletions.

Two things are not derivable and are recorded in `SparseNodePartition::emitted` (`packages/treetime/src/partition/storage/sparse.rs:21`) so that every output path reports the same sequence: a `--sample-from-profile` draw, which is one realization of the posterior rather than a property of it, and a tip, whose observed ambiguity and optional imputation depend on flags. Under the defaults, tips aside, nothing is stored and every accessor derives.

`seq.composition` is likewise no longer rewritten at output time. It describes the stored parsimony sequence, which is the accounting both `combine_messages` (which subtracts one count per variable position, keyed by that position's parsimony state) and `infer_gtr` (which pairs the counts with the edges' Fitch substitutions) assume.

### Consequence for accessors

Because the stored sequence is no longer the emitted one, every reader must go through an accessor that resolves it. There are two `node_sequence` methods - `PartitionBranchOps::node_sequence` and `AugurNodeDataJsonAncestralPartition::node_sequence` (`packages/treetime/src/partition/io/augur.rs:72`) - and both must derive. Reading the field directly emits parsimony states into the node-data JSON while the reconstructed FASTA carries MAP states.

## Numerical impact

Internal-node sequences change at two classes of position: where an ancestor's most likely state deviates from parsimony and the node itself resolved the position, and where a node inherited a deletion. Per-edge `ml_subs` follow, since they are derived from the corrected sequences, and the `parent + mutations == child` invariant continues to hold. Posteriors, the marginal log-likelihood, the inferred GTR, tip sequences, and the dense backend are unchanged.

Sparse and dense now agree byte-for-byte on the 676-tip dataset at 1000, 1500, 3000, and the full 29903 columns (950 nodes each). On `data/sc2/4500` the remaining sparse/dense differences at positions 28369 and 23009 are ones the previous implementation shared, and are tracked in [M-ancestral-sparse-dense-internal-gap-placement-diverges](../issues/M-ancestral-sparse-dense-internal-gap-placement-diverges.md).

## Test coverage

`test_marginal_map_deviation_does_not_leak_into_subtree` (`packages/treetime/src/ancestral/__tests__/test_marginal_map_deviation.rs`) builds an 8-taxon tree reproducing the first defect in miniature: a single `T` tip on a zero-length branch drives its parent's argmax to `T` at `P ~ 0.99`, while a tight bundle of `C` tips below keeps every node in that clade at `C` with `P > 1 - EPS`, so all of them resolve the position out of `variable`. It asserts the deviation stays on the node that owns it and that the whole reconstruction matches dense.

`test_marginal_map_deviation_keeps_inherited_deletions` guards the second defect's invariant but does **not** reproduce it: on that topology the deleted clade's nodes carry no posterior at the position, so it passes on the pre-parsimony-chain reconstruction too. A faithful regression test needs a node that is `non_char` at a position its own upward message still reports as variable, as at the `data/sc2/4500` positions above. This is outstanding.
