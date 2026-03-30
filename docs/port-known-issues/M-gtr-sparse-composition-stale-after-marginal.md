# Sparse GTR inference mixes MAP mutations with Fitch-era compositions

After `PartitionBranchOps` promotion, `get_mutation_counts_sparse()` derives `nij` from MAP-state mutations via `edge_subs_from_graph()`, but `Ti` base counts and `root_state` still read `SparseNodePartition.seq.composition` which is populated by Fitch compression and not updated by marginal inference.

## Impact

The `Ti` (time-in-state) base accumulates `branch_length * composition[nuc]` per edge, then adjusts for each mutation. The base counts character frequencies from Fitch compression. After marginal inference, variable-site assignments can change, making the composition stale. The adjustment terms (from MAP mutations) are correct, but the base terms are not.

Practical magnitude is small: only variable sites differ between Fitch and marginal, and these are a small fraction of total positions. The `root_state` prior is similarly affected but serves as a soft prior on equilibrium frequencies where approximate values are acceptable.

## Root cause

`SparseSeqInfo.composition` is set during Fitch compression ([packages/treetime/src/commands/ancestral/fitch.rs](../../packages/treetime/src/commands/ancestral/fitch.rs)) and never updated. `reconstruct_node_sequence()` rewrites `seq.sequence` but not `seq.composition` ([packages/treetime/src/representation/partition/marginal_sparse.rs](../../packages/treetime/src/representation/partition/marginal_sparse.rs)).

## Fix

Recompute `seq.composition` from `seq.sequence` after marginal reconstruction, or derive `Ti`/`root_state` from the same MAP state that `edge_subs_from_graph()` uses.
