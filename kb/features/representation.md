# Representation and Partition System

## Partitions

- [x] PartitionFitch (Fitch parsimony, sparse storage)
- [x] PartitionMarginalDense (full probability vectors at all positions)
- [x] PartitionMarginalSparse (variable positions only)
- [x] PartitionMarginalDiscrete (discrete trait states for mugration)
- [x] Marginal pass trait (`MarginalPasses`, shared by dense, sparse, and discrete partitions)
- [x] Partition set: partitions as sets of sites of the `--alignment` files, one per file by default or one per `--partition` entry, each with its own model; one representation per run (`enum PartitionSet`, [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md))
- [x] Partition entries with site ranges `a`, `a-b` and `a-b\k` (1-based, the syntax of RAxML-NG partition files and NEXUS charsets) and a model per entry; every site of a file belongs to one entry, and each model must fit the alphabet, also without entries ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md))
- [x] Assembly of a split alignment back to file coordinates: sequences, ambiguity codes, substitutions, and indel events counted once per file ([kb/algo/partition-assembly.md](../algo/partition-assembly.md))
- [ ] Automatic codon-position partitioning from GFF annotation input (codon positions given by hand work through `--partition` site ranges with a step of 3, for example `101-1600\3`)

## Payloads

- [x] Dense state (`DenseNodeState`, edge messages, `DenseSeqDistribution`)
- [x] Sparse state (`SparseNodeObs`, `SparseNodeState`, `SparseEdgeObs`, `SparseSeqDistribution`)
- [x] Node and edge data outside the graph, in maps keyed by `GraphNodeKey` and `GraphEdgeKey` (time distributions, branch distributions, messages: [kb/decisions/graph-topology-only-value-pipeline.md](../decisions/graph-topology-only-value-pipeline.md))

## Operations

- [x] Mutation composition (edge merge: A5G + G5T = A5T with cancellation)
- [x] Edge inversion (swap ref/qry, reverse messages for rerooting)
- [x] Reroot operations on partitions (split/merge edges, update messages)

## Architecture

- [x] Dense/sparse architecture ([intentional change](../decisions/sequence-representation-dense-sparse.md))
- [x] Partition system ([intentional change](../decisions/partition-system-architecture.md))
- [x] Graph-based phylogenetic representation ([intentional change](../decisions/graph-based-phylogenetic-representation.md))
