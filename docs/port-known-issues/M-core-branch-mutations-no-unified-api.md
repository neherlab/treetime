# Branch mutations have no unified API across partition types

After marginal reconstruction, probability distributions encode the current nucleotide state at each node (dense: full `Array2<f64>` per edge, sparse: `BTreeMap<usize, VarPos>` at variable sites). The discrete mutation on a branch is the change in argmax of the parent and child probability distributions. Sparse partitions also carry a `subs: Vec<Sub>` field populated during the earlier Fitch parsimony step, but these stored substitutions become stale once marginal inference updates the probability distributions.

The optimize command addressed this by introducing `PartitionOptimizeOps::edge_subs()` ([packages/treetime/src/commands/optimize/partition_ops.rs#L14-L27](../../packages/treetime/src/commands/optimize/partition_ops.rs#L14-L27)), a trait method that derives mutations from the current probabilistic state for both dense and sparse partitions. Other commands bypass this trait and access mutations through inconsistent paths, producing stale or missing results.

## Affected consumers

### Prune command reads stale Fitch substitutions

`get_edge_num_muts()` ([packages/treetime/src/commands/prune/run.rs#L136-L152](../../packages/treetime/src/commands/prune/run.rs#L136-L152)) reads `edge_partition.subs.len()` directly from `SparseEdgePartition`. When the prune command runs after marginal reconstruction, these counts reflect Fitch parsimony output, not the current marginal state. A branch where marginal inference resolved an ambiguity differently from Fitch will have an incorrect mutation count, affecting `--prune-empty` decisions.

### GTR inference uses different mutation representations per partition type

Sparse GTR inference `get_mutation_counts_sparse()` ([packages/treetime/src/gtr/infer_gtr/sparse.rs#L23-L69](../../packages/treetime/src/gtr/infer_gtr/sparse.rs#L23-L69)) reads `partition.edges[&edge_key].subs` directly (line 57), counting integer substitutions from the stale Fitch field. Dense GTR inference `get_mutation_counts_dense()` ([packages/treetime/src/gtr/infer_gtr/dense.rs#L127-L196](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L127-L196)) computes fractional expected counts from the full joint posterior via `get_branch_mutation_matrix()`. The two paths count biologically different quantities: integer stale-Fitch mutations vs. fractional current-posterior expectations. During iterative GTR refinement, the sparse path feeds back counts that do not reflect the model's own marginal distributions.

### Newick and Nexus output emit empty mutation annotations

`NodeAncestral::nwk_comments()` ([packages/treetime/src/representation/payload/ancestral.rs#L35](../../packages/treetime/src/representation/payload/ancestral.rs#L35)) and `NodeClock::nwk_comments()` ([packages/treetime/src/commands/clock/clock_graph.rs#L42](../../packages/treetime/src/commands/clock/clock_graph.rs#L42)) both contain `let mutations: String = "".to_owned(); // TODO: fill mutations`. The Nexus writer for timetree (`NodeTimetree::nwk_comments()` at [packages/treetime/src/representation/payload/timetree.rs#L129](../../packages/treetime/src/representation/payload/timetree.rs#L129)) builds comments from scratch with only `date`, never consulting partition data. v0 annotates branches with `[&mutations="A55G,T93C",date=2003.84]`; v1 produces only `[&date="2000.68"]`.

The partition data needed to fill these annotations exists (sparse `edge.subs` or dense probability matrices), but there is no way for the payload-layer Newick writer to access partition-layer state. This is the structural gap that the `PartitionOptimizeOps` trait solves for the optimize command but no equivalent exists at the output layer.

### Shared-mutation merging reads stale Fitch substitutions

`merge_shared_mutation_branches()` ([packages/treetime/src/commands/prune/run.rs#L311-L331](../../packages/treetime/src/commands/prune/run.rs#L311-L331)) and its helper `compute_shared_subs_across_partitions()` ([packages/treetime/src/commands/prune/run.rs#L422-L437](../../packages/treetime/src/commands/prune/run.rs#L422-L437)) read `partition.edges[&edge_key].subs` directly to find shared mutations between siblings. When the prune command runs after marginal reconstruction (or when integrated into the optimize loop), the shared-mutation detection operates on stale Fitch-era data. Siblings whose marginal posteriors differ from Fitch assignments will produce incorrect shared/unique mutation partitioning, leading to wrong topology changes and wrong branch length estimates for new internal nodes.

## Root cause

The Fitch forward pass populates `SparseEdgePartition.subs` ([packages/treetime/src/representation/payload/sparse.rs#L100](../../packages/treetime/src/representation/payload/sparse.rs#L100)) ([packages/treetime/src/commands/ancestral/fitch.rs#L442](../../packages/treetime/src/commands/ancestral/fitch.rs#L442)) and never updated after marginal inference. Dense partitions have no `subs` field at all; `DenseEdgePartition` ([packages/treetime/src/representation/payload/dense.rs#L38-L44](../../packages/treetime/src/representation/payload/dense.rs#L38-L44)) stores only probability messages.

The `PartitionOptimizeOps` trait that correctly derives mutations from current state lives in the optimize command ([packages/treetime/src/commands/optimize/partition_ops.rs#L14](../../packages/treetime/src/commands/optimize/partition_ops.rs#L14)). No equivalent trait or method exists at the partition layer for general use by other commands and output writers.

## Proposed solutions

### S1: Promote edge_subs to a general partition trait

Move `edge_subs()` from `PartitionOptimizeOps` to a trait on the partition layer (e.g. `PartitionMarginal` or a new `PartitionBranchMutations` trait in `representation/partition/traits.rs`). Both `PartitionMarginalDense` and `PartitionMarginalSparse` already implement the method; the change is organizational. All consumers (prune, GTR, Newick writers) route through the trait instead of reading `edge.subs` directly.

### S2: Update SparseEdgePartition.subs after marginal inference

After each marginal pass, overwrite `edge.subs` with the current argmax-derived mutations. Consumers that read `edge.subs` directly then get current data without API changes. Downside: the field carries dual semantics (Fitch-era vs. marginal-era) depending on when it was last written, and dense partitions still lack a `subs` field entirely, so the asymmetry persists.

### S3: S1 + retire direct subs access in post-marginal code paths

Combine S1 with a lint or convention that `edge.subs` is only read in Fitch-specific code. Post-marginal consumers use the trait. This preserves Fitch-era `subs` for compression and sparse representation while ensuring all probabilistic consumers see current state.

## v0 handling

v0 does not have this split. `TreeAnc._ml_anc()` ([packages/legacy/treetime/treetime/treeanc.py](../../packages/legacy/treetime/treetime/treeanc.py)) overwrites `node.mutations` with the current ML reconstruction result after each pass. All consumers (GTR inference, output, branch length optimization) read the same `node.mutations` dict and always see current state.

## Related issues

- [Nexus output missing mutation annotations](M-timetree-nexus-missing-mutations.md): the output-side symptom of this structural gap
- The optimize command already relies on current-state branch mutations for initial branch length estimation and gap-aware effective lengths. That existing usage shows the shared API is viable, but it remains confined to optimize-specific traits.
