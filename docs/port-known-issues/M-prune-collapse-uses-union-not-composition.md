# Prune edge collapse uses set-union instead of mutation composition

`collapse_sparse_edge()` in `packages/treetime/src/commands/prune/run.rs` merges substitutions from collapsed edges using `iterator_union()`, which performs sorted deduplication. The correct operation is mutation composition using `compose_substitutions()` (already implemented in `packages/treetime/src/representation/partition/marginal_sparse.rs`).

## Problem

When collapsing an internal node B between parent A and child C, substitutions from edges A-B and B-C must be composed to produce the net A-C substitutions:

- Non-overlapping positions: keep both (correct with union)
- Same position, chain: A0G + G0T = net A0T (union keeps both as separate entries)
- Same position, cancellation: A0G + G0A = no net change (union keeps both)
- Same position, identical: A0T + A0T (convergent mutation) = union deduplicates to one copy

Set-union satisfies neither "net state change" nor "all events" semantics.

## Impact

For `--prune-empty` (zero-mutation edges), union is accidentally correct because the collapsed edge has zero substitutions. For `--prune-short` (short edges that may have mutations), union produces incorrect substitution lists on the merged edge. The incorrect list propagates to output files and affects downstream mutation counting.

## Fix

Replace `iterator_union()` with `compose_substitutions()` in `collapse_sparse_edge()`. The function already exists and handles composition, chain, and cancellation correctly.
