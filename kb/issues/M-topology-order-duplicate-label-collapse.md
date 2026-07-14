# Target topology order silently collapses duplicate labels

`compute_target_scores()` converts the requested label sequence into a `BTreeMap<&str, usize>`. A repeated label overwrites its earlier position, so the last occurrence silently wins. [`packages/treetime-graph/src/topology_order.rs#L380-L401`](../../packages/treetime-graph/src/topology_order.rs#L380-L401)

Validation checks only that every graph leaf label occurs in the map. It does not reject duplicate target labels, duplicate leaf labels, or target labels that do not identify a leaf. [`packages/treetime-graph/src/topology_order.rs#L506-L524`](../../packages/treetime-graph/src/topology_order.rs#L506-L524)

The requested order therefore does not define a one-to-one position for every leaf. Duplicates can change mean and median subtree scores without an error, while duplicate graph leaves become indistinguishable.

## Required behavior

Validate a bijection between target-order entries and graph leaves before scoring. Report duplicate target labels, duplicate leaf labels, missing leaves, and unknown target labels separately with the offending names.

## Related issues

- [M-io-sequence-name-matching-unreliable.md](M-io-sequence-name-matching-unreliable.md)

## Related tickets

- [kb/tickets/topology-order-reject-duplicate-labels.md](../tickets/topology-order-reject-duplicate-labels.md)
