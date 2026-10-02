# Dense marginal passes read only the observed leaf gaps

Every dense marginal pass classifies the gaps and unknown characters of internal nodes from the gaps and unknown characters observed at the leaves. A pass never reads leaf gaps that an earlier pass wrote. Each pass therefore gives the same result as the first, and dense agrees with sparse, which derives its gaps from one Fitch pass with the same functions.

**Type**: Behavior change in v1 (no v0 counterpart).

**Approval**: approved by the user (2026-10-02). The rejected alternative kept the gaps that repeated passes accumulate, which depends on the number of passes.

**v1 location**: `PartitionMarginalDense` holds the leaf observations in `obs_leaves` ([packages/treetime/src/partition/marginal/dense/partition.rs](../../packages/treetime/src/partition/marginal/dense/partition.rs)). The backward pass builds every leaf from its observation in `fn leaf_backward()` ([packages/treetime/src/partition/marginal/shared/pass.rs](../../packages/treetime/src/partition/marginal/shared/pass.rs)).

## Rule

- A leaf enters every backward pass with the gaps, unknown ranges and sequence of its input sequence
- An internal node takes its gap and unknown ranges from its children, through `compute_node_ranges()` and `resolve_indels_backward()` ([packages/treetime/src/seq/indel.rs](../../packages/treetime/src/seq/indel.rs))
- The forward pass still adds a parent gap to a child that is non-character at that position (`resolve_indels_forward()`), so the emitted leaf and internal states and the edge indels are unchanged. The next backward pass does not read these forward results

## Example

Tree `((L1,L2)X,L3)root`, six columns:

| Node | Sequence |
| ---- | -------- |
| L1   | `ACNNGT` |
| L2   | `ACNNGT` |
| L3   | `AC--GT` |

Columns 2-3 of X have only unknown children, so X is unknown there and reconstructs as `ACNNGT`. Sparse gives the same result.

Before this rule, the forward pass copied the root gap at columns 2-3 into X and then into L1 and L2. The second backward pass read L1 and L2 as gapped, so X became `AC--GT`. The dense result depended on the number of passes: the nucleotide ancestral command ran two passes, the amino-acid path one, and timetree and optimize many. The test `test_marginal_consistency_leaf_unknown_under_gapped_ancestor_stays_unknown` ([`packages/treetime/src/ancestral/__tests__/test_marginal_consistency.rs`](../../packages/treetime/src/ancestral/__tests__/test_marginal_consistency.rs)) checks this example, and `test_prop_marginal_dense_update_idempotent` ([`packages/treetime/src/ancestral/__tests__/test_marginal_idempotency_prop.rs`](../../packages/treetime/src/ancestral/__tests__/test_marginal_idempotency_prop.rs)) checks that a second pass returns the result of the first.

## Why

- A pass that depends on the output of the previous pass is not a function of the data. Its result changes with the pass count, which differs between commands
- The rule is the one sparse already uses, so the two representations agree on gaps. On `data/sc2/4500` and `data/rsv/a/20` the dense and sparse reconstructions now place every gap and unknown character at the same positions
- The alternative, the fixed point of repeated passes, would also change amino-acid and deep-tree nucleotide output, and it has no single-pass definition

v0 is not a reference for this case: it models the gap as a fifth state of the `nuc` alphabet, which is a different model.

## Impact

- Output changes only for dense runs (`--dense true`) in which a leaf has an unknown character under a gapped ancestor. Internal nodes then show `N` instead of `-` at these positions
- On the smoke datasets this changed the internal sequences of `ebola/20`, `ebola/100`, `lassa/L/20`, `lassa/L/50`, `mpox/clade-ii/20` and `mpox/clade-ii/100` with `--dense true`, and every changed position now equals the sparse result
