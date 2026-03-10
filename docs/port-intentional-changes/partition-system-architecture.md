# Partition system architecture

| Property    | Value                                                                                                                                    |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Intentional architectural change                                                                                                         |
| v1 Location | [`packages/treetime/src/representation/partition/`](../../packages/treetime/src/representation/partition/) (partition types and traits)  |
| v1 Location | [`packages/treetime/src/representation/payload/`](../../packages/treetime/src/representation/payload/) (per-node/edge data structs)      |
| v0 Location | `TreeAnc` (`#TreeAnc`) in [`packages/legacy/treetime/treetime/treeanc.py`](../../packages/legacy/treetime/treetime/treeanc.py)           |
| v0 Location | `ClockTree` (`#ClockTree`) in [`packages/legacy/treetime/treetime/clock_tree.py`](../../packages/legacy/treetime/treetime/clock_tree.py) |
| v0 Location | `TreeTime` (`#TreeTime`) in [`packages/legacy/treetime/treetime/treetime.py`](../../packages/legacy/treetime/treetime/treetime.py)       |
| Affects     | Code organization, type safety, extensibility, memory layout                                                                             |
| Commands    | All commands                                                                                                                             |

## What v0 does

v0 uses a monolithic class hierarchy: `TreeTime` extends `ClockTree` extends `TreeAnc` (4411 lines total across three files). `TreeAnc` alone has 55 methods spanning 1831 lines.

The `TreeAnc` class owns:

- The Biopython tree object (`self.tree`)
- A single GTR model (`self.gtr`)
- Sequence data (`self.data`, a `SequenceData` instance)
- All ancestral reconstruction methods (Fitch, marginal, joint)
- GTR inference
- Branch length optimization

`ClockTree` adds date handling, node time distributions, and molecular clock inference on top. `TreeTime` adds rerooting, outlier detection, and polytomy resolution.

Sequence and reconstruction state is stored directly on Biopython `Clade` nodes via dynamic attribute assignment. v0 monkey-patches `Clade` at module level to add properties like `cseq` (compressed sequence) and `mutations`:

```python
Clade.sequence = property(lambda x: x.tt.sequence(x, as_string=False))
Clade.cseq = property(compressed_sequence)
Clade.mutations = property(mutations)
```

Internal nodes accumulate attributes during reconstruction (`_cseq`, `marginal_profile`, `marginal_log_Lx`, `joint_Lx`, `joint_Cx`, `branch_state`, `mask`, `count`, `dist2root`, `bad_branch`). These attributes appear and disappear at different stages of computation, with no static guarantee of their presence.

The single class handles one data type at a time. There is no mechanism for operating on multiple independent data partitions (e.g., different genes, coding vs non-coding regions) in a single tree traversal.

## What v1 does

v1 separates tree structure from per-partition data. The generic `Graph<N, E, D>` stores topology and node/edge payloads. Partition structs store reconstruction state indexed by node and edge keys.

Each partition type is a standalone struct containing its own:

- `BTreeMap<GraphNodeKey, _>` for per-node data
- `BTreeMap<GraphEdgeKey, _>` for per-edge data
- Algorithm-specific configuration (GTR model, alphabet, sequence length)

The partition types form a hierarchy matching the reconstruction pipeline:

| Partition                 | Purpose                         | Node data             | Edge data             |
| ------------------------- | ------------------------------- | --------------------- | --------------------- |
| `PartitionFitch`          | Parsimony reconstruction        | `SparseNodePartition` | `SparseEdgePartition` |
| `PartitionMarginalDense`  | Full probability matrices       | `DenseNodePartition`  | `DenseEdgePartition`  |
| `PartitionMarginalSparse` | Variable positions only         | `SparseNodePartition` | `SparseEdgePartition` |
| `PartitionLikelihood`     | GTR + alphabet config (no data) | -                     | -                     |

Behavior is defined through traits, not inheritance:

- `PartitionMarginalOps` - backward/forward passes, sequence attachment, ancestral extraction
- `PartitionCompressed` - access to sparse node/edge maps, alphabet, length
- `PartitionTimetreeOps` - branch length optimization contributions, topology reconciliation
- `PartitionRerootOps` - mutation/message updates during rerooting (no-op for dense, active for sparse)
- `HasLogLh` - log-likelihood extraction per node

The combined trait `PartitionTimetreeAll` uses a blanket implementation to compose these:

```rust
pub trait PartitionTimetreeAll<N, E>:
  PartitionMarginalOps<N, E> + PartitionTimetreeOps<N, E> + PartitionRerootOps + HasLogLh
```

Commands operate on `Vec<Arc<RwLock<dyn PartitionTimetreeAll<...>>>>`, iterating over partitions during each tree traversal. The graph log-likelihood is the sum across all partitions:

```rust
pub fn graph_log_lh<P, N, E, D>(graph: &Graph<N, E, D>, partitions: &[Arc<RwLock<P>>]) -> Result<f64, Report>
```

## Why v1 changes this

Separation of concerns. v0's `TreeAnc` mixes tree topology management, sequence data storage, multiple reconstruction algorithms, GTR inference, and branch length optimization in a single class. Each concern cannot be tested, replaced, or extended independently.

Type safety. v0 relies on dynamic attributes (`hasattr` checks, `__delattr__`, monkey-patched properties). Whether `node.marginal_profile` exists depends on which method ran last. v1's partition structs make all data statically typed - `DenseNodePartition` always has `seq: DenseSeqInfo` and `profile: DenseSeqDis`.

Multi-partition support. v0 processes one alignment at a time. v1 can attach multiple partitions to the same tree, each with its own alphabet, GTR model, and reconstruction state. Tree traversals visit each node once and update all partitions, enabling joint analysis of heterogeneous data (e.g., different genes with different substitution models).

Dense/sparse interchangeability. `PartitionMarginalDense` and `PartitionMarginalSparse` implement the same `PartitionMarginalOps` trait. Commands that accept `dyn PartitionMarginalOps` work with either representation without conditional logic.

Memory layout. v0 scatters reconstruction data across Biopython node objects as Python attributes. v1 stores partition data in contiguous `BTreeMap` collections, separate from the graph nodes, improving cache locality during tree traversals that access reconstruction state.

## Practical impact

- No user-facing behavioral difference. The partition system is an internal architectural change.
- Multiple partitions per tree enable future multi-gene analysis without re-traversing the tree per gene.
- Adding a new reconstruction method requires implementing the relevant trait(s) on a new struct, without modifying existing partition types or the graph structure.
- Dense and sparse modes are selected at partition construction time and dispatched through trait objects, so the choice is transparent to command implementations.
