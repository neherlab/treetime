# A node without a marginal state gets an empty ancestral sequence

The dense and the sparse `extract_ancestral_sequence()` return an empty sequence when a node has no marginal state:

- dense: [packages/treetime/src/partition/marginal/dense/partition.rs#L165-L175](../../packages/treetime/src/partition/marginal/dense/partition.rs#L165-L175)
- sparse: [packages/treetime/src/partition/marginal/sparse/partition.rs#L107-L116](../../packages/treetime/src/partition/marginal/sparse/partition.rs#L107-L116)

A missing node state is a broken invariant. The empty sequence reaches the output as if the node had no sites.

## Fix

Return an internal error that names the node, as the project rule for missing node data requires (crash or fail on a missing key instead of a default).

## Validation

A test removes the marginal state of one node and expects the error instead of an empty sequence.
