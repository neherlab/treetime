# Sparse branch optimization reads stale Fitch-era states

## Summary

`get_coefficients` in sparse branch-length optimization reads parent/child states from Fitch-era `msg.variable.state` fields and `fitch_subs().reff()/qry()`. After marginal inference changes MAP assignments, these states are stale and optimization scores against wrong nucleotides.

## Details

`packages/treetime/src/partition/optimize_sparse.rs`:

The sparse branch optimization builds `SiteContribution` entries using:

- `edge.msg_to_child.variable[pos].state` for parent state
- `edge.msg_to_parent.variable[pos].state` for child state
- `edge.fitch_subs()` for substitution positions

These `.state` fields are set during the Fitch forward pass and the initial `process_node_backward` call. After marginal inference (backward + forward passes), the MAP state at variable positions can change, but the `.state` fields are not updated. The optimization then computes branch-length contributions using the Fitch-era nucleotide assignments instead of the current posterior-derived assignments.

## Impact

The bias is proportional to the number of variable sites where the marginal MAP state differs from the Fitch assignment. For low-diversity viral datasets this is a small fraction, but grows with sequence divergence and ambiguity.

## Fix

Read states from the current marginal profiles (argmax of the posterior) instead of the stored `.state` fields. Alternatively, update `.state` fields after each marginal pass.
