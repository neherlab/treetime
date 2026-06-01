# Sparse marginal msg_to_child includes near-deterministic variable sites

## Summary

`process_node_forward` unconditionally inserts all variable positions into `msg_to_child.variable`, even when the MAP state matches the fixed-row fallback and the probability is near 1.0. This causes `fixed_counts` to be decremented for positions that contribute no quantitative information beyond the fixed row.

## Details

`packages/treetime/src/partition/marginal_passes.rs:387-411`:

Every variable position is inserted into `msg_to_child.variable` and its character state is decremented from `fixed_counts` (line 411). For near-deterministic positions where `max_prob >= 1.0 - EPS` and the MAP state matches the canonical parent state, the variable row is identical to the fixed row. Moving these from fixed to variable changes their treatment in branch-length optimization: the fixed path uses `fixed_counts` as a bulk multiplicity, while the variable path processes each position individually.

## Impact

The effect is small for typical data where most variable sites have genuine uncertainty. It becomes noticeable with long branches or high-confidence reconstructions where many variable positions are near-deterministic.

## Fix

Filter variable positions before insertion: only keep in the variable map when the MAP state differs from the canonical fallback or the maximum probability is below `1.0 - EPS`.
