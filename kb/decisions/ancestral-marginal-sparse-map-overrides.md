# Sparse marginal reconstruction pins states the parent chain cannot supply

This document records a correctness fix in the v1 sparse marginal ancestral reconstruction. It affects internal-node sequences emitted by `treetime ancestral --method-anc=marginal` and `treetime timetree` under the default sparse backend. Posteriors, log-likelihoods, and the dense backend are unchanged.

This is the internal-node counterpart of the tip corruption recorded in [Marginal Tip Reconstruction and Missing-Data Imputation](ancestral-marginal-tip-reconstruction-and-imputation.md). Both come from the same source - reconstructing a node by chaining off its parent - but the tip fix removed leaves from the chain, leaving the internal-node case in place.

## Background

The sparse backend does not store a posterior at every position. `combine_messages` (`packages/treetime/src/partition/marginal/sparse/message.rs:17`) keeps an explicit distribution only for positions that stay genuinely uncertain; once a position's posterior concentrates on one state above `1 - EPS` (`EPS = 1e-4`, `message.rs:15`) and every incoming message agrees with the parsimony reference state, the position is dropped from `profile.variable` and represented by the shared per-character `fixed` vector instead.

`reconstruct_map_seq_sampled` (`packages/treetime/src/partition/marginal/sparse/reconstruct.rs:17`) then rebuilds a node's sequence as: clone the parent's reconstructed sequence, apply the **parsimony** substitutions and indels recorded on the edge, and overwrite the positions the node still holds in `profile.variable`. Positions not in `variable` are, by construction, taken from the parent.

That chain carries an unstated assumption: the parent's own state at such a position equals its parsimony reference state. Where it does, the parent hands down the reference chain and the fallback is exactly right.

## Problem

The assumption breaks wherever the marginal argmax at a node picks a state Fitch did not. Marginal reconstruction is free to disagree with parsimony, and does so systematically around near-zero branch lengths: a tip on a zero-length branch is nearly conclusive evidence about its parent, so the maximum-likelihood explanation is that the parent already carried the tip's state, while parsimony places the substitution on the tip branch.

At such a node the reconstruction correctly writes the deviating state, because the position is still in its `profile.variable`. But every descendant that resolved the position out of `variable` inherits that state instead of its own. One node's legitimate state change silently becomes a subtree-wide one, and the deeper and more conserved the clade, the more nodes it flips - resolving a position out of `variable` is exactly what a deep, conserved clade does.

### Worked example

On a 676-tip SARS-CoV-2 subtree with a 1000-column alignment starting at genome position 22000, position 928 holds `C` in 674 tips and `T` in two (`JN.1.11`, `JN.1.11.2`). Tip `JN.1.11` hangs off `NODE_0000080` on a zero-length branch. Fitch places the single `C928T` on that tip branch, but the near-zero branch makes `T` at the parent more likely than a mutation: the true posterior at `NODE_0000080` is `C = 0.4680 / T = 0.5320`.

The sparse backend computed that posterior exactly - it was verified against an independent Felsenstein implementation and against the dense backend, and all three agree to every printed digit. `NODE_0000080` and `NODE_0000081` are genuinely `T`. The clade below them is all-`C` with upward evidence of order `1e-32`, so its nodes resolved position 928 out of `variable` and inherited `T` from above. The run emitted 194 internal `T` where only 2 are correct.

### Why it looked length-dependent

The apparent dependence on alignment length is indirect and was a red herring. The default `--model infer` estimates the GTR from the data, so adding alignment columns shifts the rate matrix, which moves a near-tied posterior across the 50/50 line. On this dataset the flip is caused by a single column: truncating the alignment to 985 columns reconstructs position 928 correctly, and adding column 986 (itself a `C`/`T` polymorphism) flips `NODE_0000080` and with it the whole clade. `--model jc69` never flips it. Longer alignments simply make the flip more likely and can move it higher in the tree, where the corrupted clade is larger.

Note also that zero-length branches are clamped to `MIN_BRANCH_LENGTH_FRACTION / seq_length` (`packages/treetime/src/constants.rs:7`, `packages/treetime/src/hacks/fix_branch_length.rs`), so the evidence a zero-branch tip contributes about its parent grows with alignment length. This affects both backends identically and is not part of this fix, but it is what makes these near-tied parent states common enough to matter.

Parsimony reconstruction and `--dense=true` were never affected: neither chains a node's sequence off its parent this way.

## Decision

A node records the states it cannot obtain from the parent chain, and the reconstruction applies them.

`SparseNodePartition::map_overrides` (`packages/treetime/src/partition/storage/sparse.rs:32`) is a `BTreeMap<usize, AsciiChar>` populated by the forward pass through `collect_map_overrides` (`packages/treetime/src/partition/marginal/sparse/forward.rs:162`). It holds every position that the parent still carries in `profile.variable` and that this node resolved away, mapped to the reference state the position collapsed onto. `reconstruct_map_seq_sampled` writes these with the reference chain, before the unknown and deletion masks, so gap and `N` handling still takes precedence (`reconstruct.rs:45`).

Three properties make this cheap and safe:

- **Bounded size.** Only positions the parent still holds a distribution for are candidates, so the map is bounded by the parent's `variable` set and stores one byte per entry rather than a K-wide vector. The sparse memory characteristics are unaffected.
- **One level of propagation suffices.** `compute_msg_to_child` seeds `msg_to_child.variable` from every entry of the parent's `profile.variable`, so a deviating position always reaches the immediate children as a candidate. A child that pins it back to the reference state carries no `variable` entry for it, so its own children never see the deviation.
- **No prediction of the parent's state.** The override is written unconditionally rather than only when the parent is known to deviate. It is a no-op whenever the parent agrees with parsimony. This avoids having to reproduce the parent's emitted character, which is the marginal argmax under `--sample-from-profile=argmax` but a random draw under `root` or `all` - the sampling modes are covered by the same mechanism for free.

## Numerical impact

Internal-node sequences change only at positions where an ancestor's most-likely state deviates from parsimony and the node itself resolved the position; those positions now carry the node's own state. Per-edge `ml_subs` follow, since they are derived from the corrected sequences, and the `parent + mutations == child` invariant continues to hold. Posteriors, the marginal log-likelihood, the inferred GTR, tip sequences, and the dense backend are unchanged.

Sparse and dense now agree byte-for-byte on the dataset above at 1000, 1500, 3000, and the full 29903 columns (950 nodes each, zero differing characters). Before the fix, the 1000-column run differed from dense at 204 nodes across 2 positions.

## Test coverage

`test_marginal_map_deviation_does_not_leak_into_subtree` (`packages/treetime/src/ancestral/__tests__/test_marginal_map_deviation.rs`) builds an 8-taxon tree reproducing the shape in miniature: a single `T` tip on a zero-length branch drives its parent's argmax to `T` at `P ~ 0.99`, while a tight bundle of `C` tips below keeps every node in that clade at `C` with `P > 1 - EPS`, so all of them resolve the position out of `variable`. The test asserts that the deviation stays on the one node that owns it, that the clade below stays `C`, and that the whole reconstruction matches the dense backend.
