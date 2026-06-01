# Sparse node_reference_state ignores non-char positions

## Summary

`node_reference_state()` reads `fitch.chosen_state` for variable positions but has no entry for non-variable positions with gaps or N. When a variable position overlaps a node's gap/unknown range, the returned state is a concrete nucleotide instead of gap/N.

## Details

`packages/treetime/src/partition/marginal_passes.rs`:

`node_reference_state()` returns `fitch.chosen_state.get(&pos)`. `chosen_state` is populated during the Fitch forward pass only for variable positions. At a position where the node has N or gap, `chosen_state` may contain the resolved nucleotide (Fitch resolves ambiguity to a concrete state) rather than N/gap.

`node_reference_state_or()` filters non-canonical states and falls back to `p.get_one()` from the Fitch state set. This fallback also returns a concrete nucleotide.

The `state` field in `VarPos` entries is used as the fixed-row lookup key in `msg_to_child`/`msg_to_parent`. A wrong key selects the wrong fixed-row profile, producing incorrect branch-length contributions for that position.

## Impact

Affects variable positions where a node has a gap or N. This is a small fraction of positions for typical viral datasets, but the fraction grows with divergence and ambiguity.

## Fix

Check `non_char` ranges before returning from `node_reference_state`. If the position falls in a gap range, return gap. If in an unknown range, return unknown. This matches the precedence order: gap > unknown > chosen_state > fallback.

## v0 comparison

v0 overwrites `node.cseq` after each marginal pass, embedding gap and N positions directly in the sequence. State lookup reads from `cseq`, which always reflects the correct character including gaps and N.
