# Sparse timetree convergence tracking compares empty sequences

## Problem

`capture_ancestral_states()` snapshots ancestral sequences for convergence tracking via `extract_ancestral_sequence()`. After `compress_sequences()` finalizes sparse partitions (clearing internal node sequences and storing the root sequence off-tree), `extract_ancestral_sequence()` returns empty `Seq` for non-root internal nodes. Both pre- and post-iteration snapshots are empty, so `count_sequence_changes()` returns 0, and the refinement loop may declare convergence prematurely.

## Affected paths

- `packages/treetime/src/commands/timetree/refinement.rs`: `capture_ancestral_states` at lines 51 and 92
- `packages/treetime/src/commands/timetree/convergence/sequence_changes.rs`: `capture_ancestral_states` reads `extract_ancestral_sequence` per internal node
- `packages/treetime/src/partition/marginal_sparse.rs`: `extract_ancestral_sequence` returns empty when `seq.sequence` is cleared by `finalize_fitch`

## Context

Before the off-tree root sequence refactoring, internal node sequences were retained from Fitch compression. Convergence tracking compared these Fitch-era sequences across iterations. Since marginal passes do not update `seq.sequence`, the compared sequences were already stale (reflecting Fitch resolution, not current marginal posteriors). The refactoring changed the symptom from "stale but nonzero diffs" to "empty = zero diffs".

## Fix direction

`capture_ancestral_states` should reconstruct sequences from current marginal profiles (via `ancestral_reconstruction_marginal` or `node_state_at`) instead of reading stored `seq.sequence`. This makes convergence tracking reflect actual marginal state changes between iterations.

## Related issues

- Source: [M-timetree-sparse-convergence-empty-sequences.md](../issues/M-timetree-sparse-convergence-empty-sequences.md) -- delete after full resolution
