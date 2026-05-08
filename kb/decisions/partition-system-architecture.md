# Partition system architecture

This document describes an intentional architectural change affecting code organization, type safety, extensibility, and memory layout across all commands. v1 replaces v0's monolithic class hierarchy with a partition-based system that separates tree topology from per-partition reconstruction state.

v0 source files:

- `TreeAnc` (`#TreeAnc`) in `packages/legacy/treetime/treetime/treeanc.py:50-1831:1:`
- `ClockTree` (`#ClockTree`) in `packages/legacy/treetime/treetime/clock_tree.py:11-1234:1:`
- `TreeTime` (`#TreeTime`) in `packages/legacy/treetime/treetime/treetime.py:35-800:1:`

v1 source files:

- Partition types and traits in `packages/treetime/src/representation/partition/`
- Per-node/edge data structs in `packages/treetime/src/representation/payload/`

## Background: partition models in phylogenetics

Phylogenetic trees represent evolutionary relationships inferred from molecular sequences. Ancestral sequence reconstruction estimates the sequences at internal nodes given observed sequences at leaf nodes. Two main algorithmic approaches exist: Fitch parsimony (minimizes character changes) and marginal maximum likelihood (computes probability distributions over ancestral states using a substitution model).

Substitution models describe how nucleotides evolve over time via a rate matrix Q and equilibrium frequencies. These models form a nested hierarchy from simplest to most general:

- JC69 (Jukes & Cantor, 1969): Equal base frequencies (0.25 each), equal substitution rates. Single parameter.
- K80/K2P (Kimura, 1980): Distinguishes transitions (A-G, C-T) from transversions. Equal base frequencies. Two parameters.
- HKY85 (Hasegawa, Kishino & Yano, 1985): Variable base frequencies plus transition/transversion distinction. Five parameters.
- GTR (Tavaré, 1986): General Time Reversible model with 6 exchangeability rates and 4 base frequencies. Nine free parameters. Most general time-reversible nucleotide model.

Different genomic regions (genes, codon positions) can have different evolutionary patterns. Partition models address this by allowing each region to have its own substitution model parameters while sharing the same tree topology (Lanfear et al., 2012). This is standard practice in widely-used phylogenetic software including RAxML, IQ-TREE, and BEAST.

The molecular clock hypothesis (Zuckerkandl & Pauling, 1962) proposes that mutations accumulate at roughly constant rates, enabling divergence time estimation. Treetime combines ancestral reconstruction with molecular clock inference to estimate when internal nodes existed in time.

### References

- Jukes TH, Cantor CR (1969). Evolution of protein molecules. Mammalian Protein Metabolism 3:21-132.
- Kimura M (1980). A simple method for estimating evolutionary rates of base substitutions. J Mol Evol 16:111-120.
- Hasegawa M, Kishino H, Yano T (1985). Dating of the human-ape splitting by a molecular clock of mitochondrial DNA. J Mol Evol 22:160-174.
- Tavaré S (1986). Some probabilistic and statistical problems in the analysis of DNA sequences. Lectures on Mathematics in the Life Sciences 17:57-86.
- Lanfear R, Calcott B, Ho SYW, Guindon S (2012). PartitionFinder: Combined selection of partitioning schemes and substitution models for phylogenetic analyses. Mol Biol Evol 29(6):1695-1701.
- Zuckerkandl E, Pauling L (1962). Molecular disease, evolution, and genetic heterogeneity. Horizons in Biochemistry, Academic Press, pp. 189-225.

## What v0 does

v0 uses a monolithic class hierarchy where `TreeTime` extends `ClockTree` extends `TreeAnc` (4411 lines total across three files). `TreeAnc` in `packages/legacy/treetime/treetime/treeanc.py:50-1831:1:` owns:

- The Biopython tree object (`self.tree`, assigned in `treeanc.py:165:1:`)
- A single GTR model (`self.gtr`, assigned in `treeanc.py:171:1:`)
- Sequence data (`self.data`, a `SequenceData` instance, assigned in `treeanc.py:175-184:1:`)
- All ancestral reconstruction methods: `_fitch_anc()` (`treeanc.py:575-637:1:`), `_ml_anc_marginal()` (`treeanc.py:762-812:1:`), `_ml_anc_joint()` (`treeanc.py:934-1080:1:`)
- GTR inference (`infer_gtr()`)
- Branch length optimization (`optimize_branch_lengths_joint()` in `treeanc.py:1176-1242:1:`)

`ClockTree` in `packages/legacy/treetime/treetime/clock_tree.py:11-1234:1:` adds:

- Date handling (`_assign_dates()` in `clock_tree.py:109-166:1:`)
- Node time distributions (`date_constraint`, `branch_length_interpolator`)
- Molecular clock inference (`make_time_tree()` in `clock_tree.py:401-426:1:`)
- Joint and marginal time tree optimization (`_ml_t_joint()` in `clock_tree.py:428-619:1:`, `_ml_t_marginal()` in `clock_tree.py:641-952:1:`)

`TreeTime` in `packages/legacy/treetime/treetime/treetime.py:35-800:1:` adds:

- Rerooting (`reroot()`)
- Outlier detection (`clock_filter()`)
- Polytomy resolution (`resolve_polytomies()`)
- Relaxed clock models (`relaxed_clock()`)
- The main `run()` orchestration in `treetime.py:83-450:1:`

### Dynamic attribute storage

v0 stores reconstruction state directly on Biopython `Clade` nodes via dynamic attribute assignment. Module-level monkey-patching in `treeanc.py:45-47:1:` adds properties to `Clade`:

```python
Clade.sequence = property(lambda x: x.tt.sequence(x, as_string=False))
Clade.cseq = property(compressed_sequence)
Clade.mutations = property(mutations)
```

Internal nodes accumulate attributes during reconstruction that appear and disappear at different computation stages:

- `_cseq` (compressed sequence) - assigned in `_fitch_anc()` at `treeanc.py:615-617:1:` and `_ml_anc_marginal()` at `treeanc.py:838:1:`
- `marginal_profile` - computed in `preorder_traversal_marginal()` at `treeanc.py:910-917:1:`
- `marginal_log_Lx`, `marginal_subtree_LH_prefactor` - set during `postorder_traversal_marginal()` at `treeanc.py:863-878:1:`, deleted in cleanup at `treeanc.py:804-809:1:`
- `joint_Lx`, `joint_Cx`, `seq_idx` - set during `_ml_anc_joint()` at `treeanc.py:995-996:1:`, deleted in cleanup at `treeanc.py:1069-1077:1:`
- `branch_state` - set by `add_branch_state()` at `treeanc.py:1163:1:`, conditionally deleted at `treeanc.py:891:1:` and `treeanc.py:962:1:`
- `mask`, `count`, `dist2root`, `bad_branch` - set during tree preparation at `treeanc.py:484-492:1:`

For clock tree inference, `ClockTree` adds temporal attributes:

- `raw_date_constraint`, `date_constraint` - assigned in `_assign_dates()` at `clock_tree.py:128-144:1:` and `init_date_constraints()` at `clock_tree.py:382-399:1:`
- `branch_length_interpolator` - created in `init_date_constraints()` at `clock_tree.py:360-370:1:`
- `time_before_present`, `clock_length` - assigned during `_ml_t_joint()` at `clock_tree.py:533-536:1:` and `clock_tree.py:614-615:1:`
- `joint_pos_Lx`, `joint_pos_Cx` - created during joint time tree optimization, deleted in cleanup at `clock_tree.py:449-451:1:`
- `marginal_pos_Lx`, `marginal_pos_LH`, `subtree_distribution`, `msg_from_parent` - created during marginal time tree optimization, deleted in cleanup at `clock_tree.py:669-677:1:`
- `numdate`, `date` - assigned by `convert_dates()` at `clock_tree.py:988-989:1:`

Whether an attribute exists depends on which methods have run, with no static guarantee of presence. The code uses `hasattr()` checks (e.g., `clock_tree.py:348:1:`, `clock_tree.py:350:1:`, `clock_tree.py:1003:1:`) and `__delattr__` (e.g., `treeanc.py:353:1:`) to manage this dynamic state.

### Single data type limitation

The class handles one alignment at a time. There is no mechanism for operating on multiple independent data partitions (e.g., different genes with different substitution models) in a single tree traversal.

## What v1 does

v1 separates tree structure from per-partition reconstruction state. The generic `Graph<N, E, D>` in `packages/treetime-graph/src/graph.rs` stores topology and node/edge payloads. Partition structs store reconstruction state indexed by node and edge keys, with each partition containing:

- `BTreeMap<GraphNodeKey, _>` for per-node data
- `BTreeMap<GraphEdgeKey, _>` for per-edge data
- Algorithm-specific configuration (GTR model, alphabet, sequence length)

### Partition types

The partition types in `packages/treetime/src/representation/partition/` form a hierarchy matching the reconstruction pipeline:

`PartitionFitch` in `packages/treetime/src/representation/partition/fitch.rs:8-15:1:` stores parsimony reconstruction state:

```rust
pub struct PartitionFitch {
  pub index: usize,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, SparseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, SparseEdgePartition>,
}
```

`PartitionMarginalDense` in `packages/treetime/src/representation/partition/marginal_dense.rs:23-31:1:` stores full probability matrices for all positions:

```rust
pub struct PartitionMarginalDense {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, DenseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, DenseEdgePartition>,
}
```

`PartitionMarginalSparse` in `packages/treetime/src/representation/partition/marginal_sparse.rs:24-32:1:` stores only variable positions:

```rust
pub struct PartitionMarginalSparse {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, SparseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, SparseEdgePartition>,
}
```

### Per-node and per-edge payload structures

Dense payloads in `packages/treetime/src/representation/payload/dense.rs`:

- `DenseNodePartition` (`dense.rs:16-19:1:`) contains `seq: DenseSeqInfo` (gap positions and full sequence) and `profile: DenseSeqDis` (probability matrix over all positions)
- `DenseEdgePartition` (`dense.rs:37-44:1:`) contains `msg_to_child`, `msg_to_parent`, `msg_from_child` (probability matrices for message passing) and `indels`
- `DenseSeqDis` (`dense.rs:46-52:1:`) contains `dis: Array2<f64>` (positions x states) and `log_lh: f64`

Sparse payloads in `packages/treetime/src/representation/payload/sparse.rs`:

- `SparseNodePartition` (`sparse.rs:16-19:1:`) contains `seq: SparseSeqInfo` (gap/unknown ranges, composition, Fitch state) and `profile: MarginalSparseSeqDistribution`
- `SparseEdgePartition` (`sparse.rs:98-106:1:`) contains `subs: Vec<Sub>` (substitutions), `indels: Vec<InDel>`, and directional messages
- `MarginalSparseSeqDistribution` (`sparse.rs:108-122:1:`) contains `variable: BTreeMap<usize, VarPos>` (probability vectors only at variable positions), `fixed: BTreeMap<AsciiChar, Array1<f64>>`, and `log_lh: f64`

### Trait-based behavior

Behavior is defined through traits in `packages/treetime/src/representation/partition/traits.rs` and `packages/treetime/src/commands/timetree/partition_ops.rs`, not inheritance:

`PartitionMarginalOps<N, E>` (`traits.rs:18-40:1:`) defines:

- `attach_sequences()` - initialize partition from alignment
- `process_node_backward()` - compute subtree likelihood (postorder pass)
- `process_node_forward()` - compute full marginal profile (preorder pass)
- `extract_ancestral_sequence()` - get reconstructed sequence from node
- `reconstruct_node_sequence()` - assign sequence during tree traversal

`PartitionCompressed` (`traits.rs:42-72:1:`) provides access to sparse node/edge maps, alphabet, and length.

`HasLogLh` (`traits.rs:74-77:1:`) extracts log-likelihood per node.

`PartitionTimetreeOps<N, E>` (`partition_ops.rs:29-43:1:`) defines:

- `create_edge_contribution()` - branch length likelihood for optimization
- `reconcile_topology()` - add entries for new nodes/edges after topology changes

`PartitionRerootOps` (`partition_ops.rs:14-24:1:`) handles mutation updates during rerooting. Default is no-op (dense); sparse partitions override to invert edge mutations.

The combined trait `PartitionTimetreeAll` (`partition_ops.rs:48-54:1:`) uses a blanket implementation to compose these:

```rust
pub trait PartitionTimetreeAll<N, E>:
  PartitionMarginalOps<N, E> + PartitionTimetreeOps<N, E> + PartitionRerootOps + HasLogLh
```

### Multi-partition operation

Commands operate on `Vec<Arc<RwLock<dyn PartitionTimetreeAll<...>>>>`, iterating over partitions during each tree traversal. The graph log-likelihood sums across all partitions (`traits.rs:80-96:1:`):

```rust
pub fn graph_log_lh<P, N, E, D>(graph: &Graph<N, E, D>, partitions: &[Arc<RwLock<P>>]) -> Result<f64, Report>
where
  P: HasLogLh + ?Sized,
{
  let root = graph.get_exactly_one_root()?;
  let root_key = root.read_arc().key();
  let log_lh = partitions
    .iter()
    .map(|partition| partition.read_arc().get_log_lh(root_key))
    .sum();
  Ok(log_lh)
}
```

Type aliases in `packages/treetime/src/representation/partition/timetree.rs:8-9:1:` simplify common combinations:

```rust
pub type GraphTimetree = Graph<NodeTimetree, EdgeTimetree, ()>;
pub type PartitionTimetreeAllVec = Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>>;
```

## Why v1 changes this

**Separation of concerns.** v0's `TreeAnc` mixes tree topology management, sequence data storage, multiple reconstruction algorithms, GTR inference, and branch length optimization in a single 1800-line class. Each concern cannot be tested, replaced, or extended independently. v1 separates: `Graph` owns topology; partition structs own reconstruction state; traits define operations.

**Type safety.** v0 relies on dynamic attributes (`hasattr` checks at `clock_tree.py:348:1:`, `__delattr__` at `treeanc.py:353:1:`, monkey-patched properties at `treeanc.py:45-47:1:`). Whether `node.marginal_profile` exists depends on which method ran last. v1's partition structs make all data statically typed. `DenseNodePartition` always has `seq: DenseSeqInfo` and `profile: DenseSeqDis`. The compiler enforces field presence.

**Multi-partition support.** v0 processes one alignment at a time with a single GTR model. v1 can attach multiple partitions to the same tree, each with its own alphabet, GTR model, and reconstruction state. Tree traversals visit each node once and update all partitions, enabling joint analysis of heterogeneous data (e.g., different genes with different substitution models).

**Dense/sparse interchangeability.** `PartitionMarginalDense` and `PartitionMarginalSparse` implement the same `PartitionMarginalOps` trait. Commands that accept `dyn PartitionMarginalOps` work with either representation without conditional logic. Dense stores full probability matrices; sparse stores only variable positions. The choice is made at partition construction time.

**Memory layout.** v0 scatters reconstruction data across Biopython node objects as Python attributes with arbitrary lifetimes. v1 stores partition data in contiguous `BTreeMap` collections, separate from graph nodes, improving cache locality during tree traversals that access reconstruction state.

## Practical impact

- No user-facing behavioral difference. The partition system is an internal architectural change.
- Multiple partitions per tree enable future multi-gene analysis without re-traversing the tree per gene.
- Adding a new reconstruction method requires implementing the relevant trait(s) on a new struct, without modifying existing partition types or the graph structure.
- Dense and sparse modes are selected at partition construction time and dispatched through trait objects, transparent to command implementations.
