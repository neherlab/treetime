# Partition trait module groups unrelated capabilities

`partition/traits.rs` contains GTR ownership, graph topology access, mutation extraction and Newick comments, marginal passes, compressed storage, transition counting, optimization, rerooting, timetree composition, and graph likelihood aggregation [`packages/treetime/src/partition/traits.rs`](../../packages/treetime/src/partition/traits.rs).

The individual traits are often narrow, but the central module and `PartitionTimetreeAll` composite make unrelated capability changes co-locate in a high-afferent boundary [`packages/treetime/src/partition/traits.rs#L335`](../../packages/treetime/src/partition/traits.rs#L335).

## Coupling mechanism

Consumers import the category module for capabilities owned by different domains. The all-capabilities composite then requires one partition type to satisfy ancestral propagation, optimization, rerooting, timetree, mutation, and likelihood concerns together. Adding a method for one workflow can expand bounds and imports for every composite consumer.

## Required organization

- Place capability traits beside the domain operation that defines their semantics.
- Keep storage-specific access beside its concrete partition owner.
- Define composite aliases at the orchestration boundary that needs the exact set.
- Prefer role-specific bounds over a crate-wide `all` trait in reusable domain functions.

This is a module and dependency correction. It does not authorize changing mathematical behavior or the approved partition architecture.

## Validation

- Dependency analysis shows domain modules import only the capabilities they use.
- Each partition implementation's trait set remains explicit and compiler-checked.
- Ancestral, optimize, reroot, timetree, and output tests preserve behavior.
- No `mod.rs` gains implementations or re-exports while modules are reorganized.

## Related issues

- [M-partition-compressed-exposes-fitch-storage.md](M-partition-compressed-exposes-fitch-storage.md)
- [N-gtr-site-specific-partition-integration.md](N-gtr-site-specific-partition-integration.md)
