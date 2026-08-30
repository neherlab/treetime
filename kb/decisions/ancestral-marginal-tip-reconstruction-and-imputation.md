# Marginal Tip Reconstruction and Missing-Data Imputation

This document records a correctness fix and a restored v0 behavior in the v1 Rust marginal ancestral reconstruction, together with the flag surface that controls tip output. The changes affect how leaf (tip) sequences are emitted by `treetime ancestral --method-anc=marginal` and `treetime timetree`. Internal-node and root reconstruction are unchanged.

## Background

Marginal reconstruction runs a backward and a forward message-passing pass over the tree, producing a per-site posterior at every node. The v1 sparse backend compresses invariant sites via Fitch parsimony and stores only variable positions densely; the dense backend stores a full posterior matrix. Both backends emit a reconstructed sequence per node, read back by the node-data JSON serializer and the reconstructed-FASTA writer.

## Problem: tip corruption (C3)

The sparse backend reconstructed every non-root node, leaves included, by cloning the parent's reconstructed sequence and re-applying the edge's substitutions and indels, then overwriting only the node's own posterior-variable positions. A leaf that is Fitch-equal to its parent at a site carries no variable entry there, so it inherited the parent's most-likely state and lost its own observed nucleotide.

Concretely, on `(A:0.4,B:0.1)` with `A=ACGT` and `B=GCGT`, the root's most-likely state at position 0 is `G`. The sparse tip `A` was emitted as `GCGT` even though its reported mutation `G1A` was correct. The mutation record and the emitted sequence disagreed.

The same overwrite ran before the per-edge maximum-likelihood substitutions (`ml_subs`) were computed, so on some edges the corrupted tip also dropped the real parent-to-tip substitution from the reported mutations. This broke the invariant that applying an edge's mutations to the parent sequence reproduces the child sequence.

The dense backend never had this defect: it keeps each leaf's observed sequence and derives leaf mutations by comparing most-likely states.

## Problem: missing data never resolved (C1)

Ambiguous (`N`, IUPAC) and unknown tip positions were stored as unresolved and re-filled with the unknown character. They were never resolved to an inferred state, regardless of any flag. v0 resolves them: its marginal leaf pass computes a full posterior at each tip and assigns the most-likely state, replacing ambiguous characters with the inferred value (`def TreeAnc.preorder_traversal_marginal()` in `packages/legacy/treetime/treetime/treeanc.py`).

## Decision

### Tips reconstruct from their own state

Leaves are reconstructed from their own observed data, not by chaining through the parent. The sparse forward pass no longer overwrites a leaf's stored sequence with the parent-chained reconstruction; the leaf keeps its observed input, mirroring the dense backend. This restores the observed tip states (C3) and restores the parent-to-tip substitutions the corruption dropped, so `parent + mutations == child` holds on every edge and the two backends agree.

### Imputation restores v0 behavior

With imputation enabled, every ambiguous or unknown tip position resolves to the argmax of the tip's marginal posterior, leaving gaps as gaps (a gap is inferred structure, not missing data). The posterior is the parent marginal evolved across the tip's branch, restricted by the observed ambiguity mask, matching v0's per-leaf marginal profile. The parent state is not simply copied, because on a long tip branch the evolved posterior can differ from the parent's most-likely state.

The dense backend takes the argmax of the stored per-site leaf profile at imputed positions; the sparse backend evolves the stored parent-to-child message across the branch and takes the argmax under the observed ambiguity mask. The two backends produce identical tip sequences.

### Flag surface

Two orthogonal flags control tip output on both `ancestral` and `timetree`:

- `--include-leaves` emits reconstructed tip sequences in addition to internal nodes.
- `--impute-missing-data` resolves ambiguous and unknown tip states to the most likely inferred state.

`--reconstruct-tip-states` is retained as the v0-compatible alias enabling both. The app-server request contracts expose all three fields with the same alias semantics. Under `--method-anc=parsimony`, imputation has no posterior to draw from; `--impute-missing-data` is a no-op that logs a warning rather than failing, so a caller that always passes `--reconstruct-tip-states` still runs, and `--include-leaves` still emits observed tips.

## Consistency of output paths

The node-data JSON serializer reads each partition's stored sequence. Reconstruction now writes that stored sequence for both backends (the dense backend previously re-derived tips from the profile and re-stamped the unknown character on every serialization), so the JSON, the reconstructed FASTA, and the sparse and dense backends all agree. Tips are always reconstructed into the stored sequence so the serializer reflects the corrected, flag-aware result; `--include-leaves` gates only whether tips are emitted to the reconstructed-FASTA stream.

## Numerical impact

Emitted tip sequences change: previously corrupted canonical tip states are restored to their observed values, and, under imputation, ambiguous and unknown tip positions resolve to the marginal argmax (matching v0). Internal-node and root sequences, per-internal-edge mutations, and the marginal log-likelihood are unchanged. Tip-edge mutations change only where the previous corruption dropped a real substitution; they now satisfy `parent + mutations == child` and match the dense backend.
