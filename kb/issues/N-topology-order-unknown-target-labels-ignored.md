# Target topology order ignores labels that match no leaf

The target order of `--topology-order=target-order` and `target-order-reverse` is a position for each graph leaf, keyed by node. It comes from one of three sources ([`packages/app-commands/src/commands/shared/topology_order_args.rs#L63-L141`](../../packages/app-commands/src/commands/shared/topology_order_args.rs#L63-L141)):

- Input order: the leaf order of the input tree as node keys, without names. `prune` captures it before pruning ([`packages/app-commands/src/commands/prune/run.rs#L45`](../../packages/app-commands/src/commands/prune/run.rs#L45)) and applies it after ([`packages/app-commands/src/commands/prune/run.rs#L100-L104`](../../packages/app-commands/src/commands/prune/run.rs#L100-L104)). Keys of pruned leaves remain in the order and match no graph leaf
- `--topology-order-tree` and `--topology-order-list`: leaf names, paired with graph leaves by `fn target_positions()` through `pair_by_name()` ([`packages/app-commands/src/commands/shared/topology_order_args.rs#L143-L153`](../../packages/app-commands/src/commands/shared/topology_order_args.rs#L143-L153)). The first position of a repeated target name wins, and every leaf with that name receives it. The test `test_input_pairing_target_order_gives_leaves_that_share_a_name_its_first_position` asserts this ([`packages/app-commands/src/commands/shared/__tests__/test_input_pairing.rs#L121-L131`](../../packages/app-commands/src/commands/shared/__tests__/test_input_pairing.rs#L121-L131))

`fn validate_target_order()` rejects a graph leaf that is absent from the target order and names it ([`packages/treetime-graph/src/topology_order.rs#L369-L385`](../../packages/treetime-graph/src/topology_order.rs#L369-L385)). It never checks the reverse. `fn target_positions()` keeps only `NamePairing::by_node` and drops `unmatched` (target names that match no leaf) and `duplicate_entry_names` (repeated target names) without a warning. The test `topology_order_target_order_ignores_absent_ranking_labels` asserts that an absent label is ignored, with a target order that contains the label `removed` ([`packages/treetime-graph/src/__tests__/test_topology_order.rs#L362-L375`](../../packages/treetime-graph/src/__tests__/test_topology_order.rs#L362-L375)).

For the input order, unmatched keys are expected, because they belong to pruned leaves. For the reference-topology and list sources, an unmatched name is a typo, a foreign label, or a leaf that the command removed.

> [!IMPORTANT]
> **Decision required.** What should happen to target names of `--topology-order-tree` and `--topology-order-list` that match no graph leaf (`unmatched`) or occur more than once (`duplicate_entry_names`)? Options: ignore them (current behavior), raise a `RunWarning` that names them, as `warn_duplicate_names()` does for other inputs, or reject them with an error. Rejecting unmatched names breaks a list or reference file that names leaves the command prunes, and contradicts `topology_order_target_order_ignores_absent_ranking_labels`. Ignoring them lets a misspelled or foreign label pass silently. Evidence: [`packages/app-commands/src/commands/shared/topology_order_args.rs#L143-L153`](../../packages/app-commands/src/commands/shared/topology_order_args.rs#L143-L153), [`packages/treetime-graph/src/pair_by_name.rs#L8-L49`](../../packages/treetime-graph/src/pair_by_name.rs#L8-L49).

## Test coverage gap

- No test reads a target order from a list file or a reference topology with a name that matches no leaf, or with a repeated name, through the command

## Related issues

- [N-io-sequence-name-matching-unreliable.md](N-io-sequence-name-matching-unreliable.md)
