# Representation and Partition System

## Partitions

- [x] PartitionFitch (Fitch parsimony, sparse storage)
- [x] PartitionMarginalDense (full probability vectors at all positions)
- [x] PartitionMarginalSparse (variable positions only)
- [x] PartitionMarginalDiscrete (discrete trait states for mugration)
- [x] Marginal pass trait (`MarginalPasses`, shared by dense, sparse, and discrete partitions)
- [x] Partition set: one partition per `--alignment` file, one representation per run (`enum PartitionSet`, [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md))
- [ ] Codon-position partitioning (1st/2nd/3rd position splits, requires GFF annotation input, probabilistic methods only)

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
