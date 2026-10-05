# Target topology order ignores labels that match no leaf

Target-order validation now rejects duplicate target labels ([`packages/treetime-graph/src/topology_order.rs#L293-L300`](../../packages/treetime-graph/src/topology_order.rs#L293-L300)), graph leaves absent from the target order ([`packages/treetime-graph/src/topology_order.rs#L441-L443`](../../packages/treetime-graph/src/topology_order.rs#L441-L443)), and duplicate graph leaf labels ([`packages/treetime-graph/src/topology_order.rs#L445-L450`](../../packages/treetime-graph/src/topology_order.rs#L445-L450)). Each error names the offending label.

Target labels that do not identify a graph leaf are accepted and ignored. `validate_target_order()` checks only that each leaf occurs in the target order, never the reverse ([`packages/treetime-graph/src/topology_order.rs#L431-L452`](../../packages/treetime-graph/src/topology_order.rs#L431-L452)). The test `topology_order_target_order_ignores_absent_ranking_labels` asserts this behavior with a target order that contains the label `removed` ([`packages/treetime-graph/src/__tests__/test_topology_order.rs#L326-L342`](../../packages/treetime-graph/src/__tests__/test_topology_order.rs#L326-L342)). The target order comes from one of three sources ([`packages/app-commands/src/commands/shared/topology_order_args.rs#L117-L146`](../../packages/app-commands/src/commands/shared/topology_order_args.rs#L117-L146)): the input leaf order, a reference topology, or a list file. The input leaf order can contain leaves that the command removes later: `prune` captures it before pruning ([`packages/app-commands/src/commands/prune/run.rs#L48`](../../packages/app-commands/src/commands/prune/run.rs#L48)) and applies it after ([`packages/app-commands/src/commands/prune/run.rs#L102`](../../packages/app-commands/src/commands/prune/run.rs#L102)).

> [!IMPORTANT]
> **Decision required.** Should target labels that match no graph leaf be rejected? A strict bijection between target-order entries and graph leaves reports every unknown label, which catches typos and wrong reference files, but it breaks the input-order source whenever a command prunes leaves, and it contradicts the existing test. Accepting unknown labels (current behavior) keeps pruned inputs working but lets a misspelled or foreign label pass silently. A per-source rule is also possible: ignore unknown labels for the input order, reject them for the reference-topology and list sources. Evidence: [`packages/treetime-graph/src/topology_order.rs#L431-L452`](../../packages/treetime-graph/src/topology_order.rs#L431-L452), [`packages/treetime-graph/src/__tests__/test_topology_order.rs#L326-L342`](../../packages/treetime-graph/src/__tests__/test_topology_order.rs#L326-L342).

## Test coverage gap

The negative tests cover one duplicate target label in the middle of the order with mean aggregation, and one duplicate graph leaf label ([`packages/treetime-graph/src/__tests__/test_topology_order.rs#L291-L324`](../../packages/treetime-graph/src/__tests__/test_topology_order.rs#L291-L324)). Missing coverage:

- Parameterized negative tests for every validation class, each asserting that the diagnostic names the offending label, including leaves absent from the target order
- Duplicate target labels at the first, middle, and last positions, under both mean and median aggregation
- Duplicates read from each target source: input order, list file, and reference topology
- Bijection validation tested independently of aggregation

## Related issues

- [M-io-sequence-name-matching-unreliable.md](M-io-sequence-name-matching-unreliable.md)
