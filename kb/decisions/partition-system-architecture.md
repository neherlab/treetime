# Partition system architecture

This document describes an intentional architectural change affecting code organization, type safety, extensibility, and memory layout across all commands. v1 replaces v0's monolithic class hierarchy with a partition-based system that separates tree topology from per-partition reconstruction state.

v0 source files:

- `TreeAnc` (`#TreeAnc`) in `packages/legacy/treetime/treetime/treeanc.py:50-1831:1:`
- `ClockTree` (`#ClockTree`) in `packages/legacy/treetime/treetime/clock_tree.py:11-1234:1:`
- `TreeTime` (`#TreeTime`) in `packages/legacy/treetime/treetime/treetime.py:35-800:1:`

v1 source files:

- Partition types, the partition set, and the marginal passes in `packages/treetime/src/partition/`
- Per-node and per-edge state structs in `packages/treetime/src/partition/storage/`

## Background: partition models in phylogenetics

Phylogenetic trees represent evolutionary relationships inferred from molecular sequences. Ancestral sequence reconstruction estimates the sequences at internal nodes given observed sequences at leaf nodes. Two main algorithmic approaches exist: Fitch parsimony (minimizes character changes) and marginal maximum likelihood (computes probability distributions over ancestral states using a substitution model).

Substitution models describe how nucleotides evolve via a rate matrix Q and equilibrium frequencies. These models form a nested hierarchy from simplest to most general:

- JC69 (Jukes & Cantor, 1969): Equal base frequencies (0.25 each), equal substitution rates. Single parameter.
- K80/K2P (Kimura, 1980): Distinguishes transitions (A-G, C-T) from transversions. Equal base frequencies. Two parameters.
- HKY85 (Hasegawa, Kishino & Yano, 1985): Variable base frequencies plus transition/transversion distinction. Five parameters.
- GTR (Tavaré, 1986): General Time Reversible model with 6 exchangeability rates and 4 base frequencies. Nine free parameters. Most general time-reversible nucleotide model.

Different genomic regions (genes, codon positions) can have different evolutionary patterns. Partition models address this by allowing each region to have its own substitution model parameters while sharing the same tree topology (Lanfear et al., 2012). This is standard practice in widely used phylogenetic software including RAxML, IQ-TREE, and BEAST.

The molecular clock hypothesis (Zuckerkandl & Pauling, 1962) proposes that mutations accumulate at roughly constant rates, enabling divergence time estimation. Treetime combines ancestral reconstruction with molecular clock inference to estimate when internal nodes existed in time.

### References

- Jukes, Thomas H., and Charles R. Cantor. 1969. "Evolution of Protein Molecules." In _Mammalian Protein Metabolism,_ vol. 3, edited by H. N. Munro, 21-132. Academic Press. https://doi.org/10.1016/B978-1-4832-3211-9.50009-7
- Kimura, Motoo. 1980. "A Simple Method for Estimating Evolutionary Rates of Base Substitutions Through Comparative Studies of Nucleotide Sequences." _Journal of Molecular Evolution_ 16(2):111-120. https://doi.org/10.1007/BF01731581
- Hasegawa, Masami, Hirohisa Kishino, and Taka-aki Yano. 1985. "Dating of the Human-Ape Splitting by a Molecular Clock of Mitochondrial DNA." _Journal of Molecular Evolution_ 22(2):160-174. https://doi.org/10.1007/BF02101694
- Tavaré, Simon. 1986. "Some Probabilistic and Statistical Problems in the Analysis of DNA Sequences." _Lectures on Mathematics in the Life Sciences_ 17:57-86.
- Lanfear, Robert, Brett Calcott, Simon Y. W. Ho, and Stephane Guindon. 2012. "PartitionFinder: Combined Selection of Partitioning Schemes and Substitution Models for Phylogenetic Analyses." _Molecular Biology and Evolution_ 29(6):1695-1701. https://doi.org/10.1093/molbev/mss020
- Zuckerkandl, Emile, and Linus Pauling. 1962. "Molecular Disease, Evolution, and Genetic Heterogeneity." In _Horizons in Biochemistry,_ edited by Michael Kasha and Bernard Pullman, 189-225. Academic Press.

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

v1 separates tree structure from per-partition reconstruction state. `struct Graph` in `packages/treetime-graph/src/graph.rs` stores the topology, branch lengths are maps keyed by `GraphEdgeKey`, and partition structs store reconstruction state in `BTreeMap` collections keyed by node and edge keys. Each partition holds:

- The observations of its sites at the nodes and edges
- Its alphabet and its number of sites
- In a reconstruction, its GTR model and the node states and edge messages of the last marginal update

### Partition types

The partition types in `packages/treetime/src/partition/` follow the reconstruction pipeline:

- `struct PartitionFitch` ([packages/treetime/src/partition/fitch/partition.rs#L17](../../packages/treetime/src/partition/fitch/partition.rs#L17)): Fitch parsimony state of one sequence partition. Its conversions produce the marginal partitions
- `struct PartitionMarginalDense` ([packages/treetime/src/partition/marginal/dense/partition.rs#L31](../../packages/treetime/src/partition/marginal/dense/partition.rs#L31)): leaf observations at all sites
- `struct PartitionMarginalSparse` ([packages/treetime/src/partition/marginal/sparse/partition.rs#L28](../../packages/treetime/src/partition/marginal/sparse/partition.rs#L28)): the root sequence and the observations at variable sites, and the total length of the run for the branch-length floor
- `struct PartitionMarginalDiscrete` ([packages/treetime/src/partition/marginal/discrete/partition.rs#L14](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L14)): leaf observations of one discrete trait (mugration)

`struct DenseReconstruction` and `struct SparseReconstruction` ([packages/treetime/src/partition/marginal/reconstruction.rs#L317-L423](../../packages/treetime/src/partition/marginal/reconstruction.rs#L317-L423)) join a sequence partition with its GTR model, its node states, and its edge messages. `enum MarginalReconstruction` holds one of the two for a single partition, and `enum ReconstructionRef` borrows one of them for the operations that both representations have, for example the log-likelihood, the edge substitutions, the edge indels, and the branch-length contributions.

The per-node and per-edge state structs are in `packages/treetime/src/partition/storage/`: `dense.rs` (`DenseNodeState`, `DenseSeqDistribution`, edge messages), `sparse.rs` (`SparseNodeObs`, `SparseNodeState`, `SparseEdgeObs`, `SparseSeqDistribution` with the profiles of the variable sites and per-character profiles of the fixed sites), and `discrete.rs`.

### Marginal passes

`trait MarginalPasses` ([packages/treetime/src/partition/marginal/shared/update.rs#L10](../../packages/treetime/src/partition/marginal/shared/update.rs#L10)) defines the backward pass, the forward pass, and the transition counts for GTR inference, with associated types for the node states and the messages of each representation. `PartitionMarginalSparse` implements it directly. The dense and discrete partitions implement it through `trait IndexedPartition` ([packages/treetime/src/partition/marginal/shared/pass.rs#L364](../../packages/treetime/src/partition/marginal/shared/pass.rs#L364)), because both store a profile matrix per node. Generic code such as GTR refinement and the rate step (`fn optimize_gtr_rate()` in `packages/treetime/src/gtr/refinement.rs`) takes any `P: MarginalPasses`, so the compiler resolves each call to one representation.

### Multi-partition operation

`enum PartitionSet` ([packages/treetime/src/partition/set.rs#L49-L52](../../packages/treetime/src/partition/set.rs#L49-L52)) holds the sequence partitions of a run in the order of the `--alignment` files:

```rust
pub enum PartitionSet {
  Dense(Vec<DenseReconstruction>),
  Sparse(Vec<SparseReconstruction>),
}
```

All partitions of a set use one representation, and construction rejects a mix. `fn build_partition_set()` ([packages/treetime/src/partition/set.rs#L19-L46](../../packages/treetime/src/partition/set.rs#L19-L46)) builds one reconstruction per file with `fn build_marginal_partition()` and gives every partition the total length of the set. The set owns the operations over all partitions: the marginal update returns the updated set and the summed log-likelihood, and the set sums the edge substitution counts, reconciles all partitions after a reroot or a topology change, and reads, sets, and estimates the relative rates. Optimize and timetree pass one set through their loops, so each loop step updates every partition with the same branch lengths. Ancestral reconstructs each file as its own `MarginalReconstruction`, because its partitions share no estimated parameters. [kb/decisions/multi-alignment-partitions.md](multi-alignment-partitions.md) records the rules for several files.

## Why v1 changes this

**Separation of concerns.** v0's `TreeAnc` mixes tree topology management, sequence data storage, multiple reconstruction algorithms, GTR inference, and branch length optimization in a single 1800-line class. Each concern cannot be tested, replaced, or extended independently. v1 separates: `Graph` owns topology; partition structs own reconstruction state; the partition set owns the operations over all partitions.

**Type safety.** v0 relies on dynamic attributes (`hasattr` checks at `clock_tree.py:348:1:`, `__delattr__` at `treeanc.py:353:1:`, monkey-patched properties at `treeanc.py:45-47:1:`). Whether `node.marginal_profile` exists depends on which method ran last. v1's partition structs make all data statically typed, and a reconstruction holds its node states and edge messages as fields, so the compiler enforces their presence.

**Multi-partition support.** v0 processes one alignment at a time with a single GTR model. v1 attaches several partitions to the same tree, each with its own alphabet, GTR model, and reconstruction state, which enables joint analysis of heterogeneous data, for example genome segments or genes with different substitution models.

**Dense/sparse interchangeability.** The dense and sparse representations implement the same `trait MarginalPasses`, and `enum PartitionSet` and `enum ReconstructionRef` dispatch the operations that commands use. Commands work with either representation without conditional logic. Dense stores full probability matrices; sparse stores only variable positions. The choice is made at partition construction time, once per run.

**Memory layout.** v0 scatters reconstruction data across Biopython node objects as Python attributes with arbitrary lifetimes. v1 stores partition data in `BTreeMap` collections separate from the graph nodes, and each marginal update returns new state instead of changing shared state in place.

## Practical impact

- The partition system is an internal architectural change. Single-file runs have no user-facing behavioral difference from it
- Several alignment files give several partitions on one tree ([kb/decisions/multi-alignment-partitions.md](multi-alignment-partitions.md))
- Adding a new reconstruction method requires a new partition struct that implements `trait MarginalPasses`, without changes to the existing partition types or the graph structure
- Dense and sparse modes are selected at partition construction time and dispatched through `enum PartitionSet` and `enum ReconstructionRef`, transparent to command implementations
