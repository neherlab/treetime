# Representation and Partition System

## Partitions

- [x] PartitionFitch (Fitch parsimony, sparse storage)
- [x] PartitionMarginalDense (full probability vectors at all positions)
- [x] PartitionMarginalSparse (variable positions only)
- [x] Partition traits (PartitionMarginal, PartitionMarginalOps, PartitionCompressed, HasLogLh)
- [ ] Codon-position partitioning (1st/2nd/3rd position splits, requires GFF annotation input, probabilistic methods only)

## Payloads

- [x] Dense payloads (DenseNodePartition, DenseEdgePartition, DenseSeqDistribution)
- [x] Sparse payloads (SparseNodePartition, SparseEdgePartition, SparseSeqDistribution)
- [x] Node payloads (NodeAncestral, NodeTimetree with time distributions)
- [x] Edge payloads (EdgeAncestral, EdgeTimetree with branch distributions and messages)

## Operations

- [x] Mutation composition (edge merge: A5G + G5T = A5T with cancellation)
- [x] Edge inversion (swap ref/qry, reverse messages for rerooting)
- [x] Reroot operations on partitions (split/merge edges, update messages)

## Architecture

- [x] Dense/sparse architecture ([intentional change](../decisions/sequence-representation-dense-sparse.md))
- [x] Partition system ([intentional change](../decisions/partition-system-architecture.md))
- [x] Graph-based phylogenetic representation ([intentional change](../decisions/graph-based-phylogenetic-representation.md))
