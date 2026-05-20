# Sparse substitution accessors

Production call sites for `ml_subs` and `fitch_subs` on `SparseEdgePartition`, excluding tests.

Definitions: [packages/treetime/src/partition/sparse.rs](../../packages/treetime/src/partition/sparse.rs)

## ML subs writes (5)

### Store computed ML subs after forward pass

After the child profile is finalized during the forward (root-to-tip) marginal pass, compute the ML substitutions for each parent edge by comparing parent and child marginal profiles, then store them on the edge.

- `process_node_forward()` [partition/marginal_passes.rs#L250](../../packages/treetime/src/partition/marginal_passes.rs#L250)
- `compute_ml_subs_for_edge()` [partition/marginal_passes.rs#L339](../../packages/treetime/src/partition/marginal_passes.rs#L339)
- `edge.set_ml_subs(ml_subs)` [partition/marginal_passes.rs#L341](../../packages/treetime/src/partition/marginal_passes.rs#L341)

### Invalidate ML subs on edge split during reroot

The old edge moves to the child side of a new internal node. ML subs are cleared because they were computed under the old root and are stale.

- `apply_reroot()` [partition/marginal_sparse.rs#L133](../../packages/treetime/src/partition/marginal_sparse.rs#L133)
- `old_edge_data.clear_ml_subs()` [partition/marginal_sparse.rs#L146](../../packages/treetime/src/partition/marginal_sparse.rs#L146)

### Invalidate ML subs after edge inversion during reroot

After inverting fitch subs on edges along the old-root-to-new-root path, ML subs are cleared because they depend on directionality that just changed.

- `apply_reroot()` [partition/marginal_sparse.rs#L133](../../packages/treetime/src/partition/marginal_sparse.rs#L133)
- `edge_data.clear_ml_subs()` [partition/marginal_sparse.rs#L212](../../packages/treetime/src/partition/marginal_sparse.rs#L212)

## ML subs reads (3)

### Count ML subs for initial branch length guess

Count discrete substitutions per edge across all partitions to estimate initial branch lengths before the optimization loop. Called from `optimize` and `timetree` commands.

- `initial_guess_mixed()` [optimize/dispatch.rs#L729](../../packages/treetime/src/optimize/dispatch.rs#L729)
- `partition.edge_subs()` [optimize/dispatch.rs#L774](../../packages/treetime/src/optimize/dispatch.rs#L774)
- `edge.ml_subs()` [partition/marginal_sparse.rs#L269](../../packages/treetime/src/partition/marginal_sparse.rs#L269)

### Collect ML subs for branch mutation annotation

Collect all ML subs across partitions for each edge and write them as comma-separated mutation strings onto the graph nodes. Used for output in `ancestral`, `optimize`, and `timetree` commands.

- `annotate_branch_mutations()` [payload/ancestral.rs#L152](../../packages/treetime/src/payload/ancestral.rs#L152)
- `partition.edge_subs()` [payload/ancestral.rs#L170](../../packages/treetime/src/payload/ancestral.rs#L170)
- `edge.ml_subs()` [partition/marginal_sparse.rs#L269](../../packages/treetime/src/partition/marginal_sparse.rs#L269)

## Fitch subs call sites (22)

### Append resolved subs during Fitch forward pass

During the forward (tip-to-root) Fitch parsimony pass, substitutions are resolved by comparing child and parent Fitch sets at each variable position. Resolved subs are appended to the edge.

- `run_fitch_forward()` [ancestral/fitch.rs#L292](../../packages/treetime/src/ancestral/fitch.rs#L292)
- `edge.extend_fitch_subs(subs)` [ancestral/fitch.rs#L453](../../packages/treetime/src/ancestral/fitch.rs#L453)

### Apply fitch subs to reconstruct child sequence for output

After Fitch reconstruction, apply fitch subs to the parent sequence to produce the child node's reconstructed sequence for output.

- `run_fitch_reconstruction()` [ancestral/fitch.rs#L561](../../packages/treetime/src/ancestral/fitch.rs#L561)
- `edge_part.fitch_subs()` [ancestral/fitch.rs#L580](../../packages/treetime/src/ancestral/fitch.rs#L580)

### Seed candidate positions for ML sub computation

Fitch subs provide candidate positions (union with parent/child variable sites) where ML substitutions may differ from Fitch. The ML computation refines these using posterior probabilities.

- `compute_ml_subs_for_edge()` [partition/marginal_passes.rs#L81](../../packages/treetime/src/partition/marginal_passes.rs#L81)
- `.fitch_subs()` [partition/marginal_passes.rs#L100](../../packages/treetime/src/partition/marginal_passes.rs#L100)

### Identify variable positions in backward marginal pass

During the backward (tip-to-root) marginal pass, fitch subs identify variable positions on each child edge. These positions seed the parent's variable-site map for marginal profile computation.

- `process_node_backward()` [partition/marginal_passes.rs#L128](../../packages/treetime/src/partition/marginal_passes.rs#L128)
- `edge_data.fitch_subs()` [partition/marginal_passes.rs#L174](../../packages/treetime/src/partition/marginal_passes.rs#L174)

### Populate state maps in forward marginal pass

During the forward (root-to-tip) marginal pass, fitch subs on parent edges identify positions where the parent and child states differ. These positions populate the parent/child state maps for marginal profile refinement.

- `process_node_forward()` [partition/marginal_passes.rs#L250](../../packages/treetime/src/partition/marginal_passes.rs#L250)
- `edge_data.fitch_subs()` [partition/marginal_passes.rs#L272](../../packages/treetime/src/partition/marginal_passes.rs#L272)

### Provide state pairs for child message precalculation

During message precalculation in the forward pass, fitch subs on each child edge provide parent/child state pairs. These are combined with sibling messages to compute the outgoing message to each child.

- `process_node_forward()` [partition/marginal_passes.rs#L250](../../packages/treetime/src/partition/marginal_passes.rs#L250)
- `child_edge_data.fitch_subs()` [partition/marginal_passes.rs#L364](../../packages/treetime/src/partition/marginal_passes.rs#L364)

### Update root sequence from parent edge during reroot merge

During edge merge in reroot, fitch subs on the parent edge update the root sequence to the old root's state (applying ref states at each sub position).

- `apply_reroot()` [partition/marginal_sparse.rs#L133](../../packages/treetime/src/partition/marginal_sparse.rs#L133)
- `edge.fitch_subs()` [partition/marginal_sparse.rs#L171](../../packages/treetime/src/partition/marginal_sparse.rs#L171)

### Compose parent and child subs during reroot edge merge

When reroot removes an internal node, the parent and child edges are merged into one. The two substitution lists are composed into a single consistent list for the merged edge.

- `apply_reroot()` [partition/marginal_sparse.rs#L133](../../packages/treetime/src/partition/marginal_sparse.rs#L133)
- `parent_edge.chain_fitch_subs(child_edge.fitch_subs())` [partition/marginal_sparse.rs#L189](../../packages/treetime/src/partition/marginal_sparse.rs#L189)

### Store composed subs on merged edge after reroot

Store the composed substitution list on the newly created merged edge after edge merge.

- `apply_reroot()` [partition/marginal_sparse.rs#L133](../../packages/treetime/src/partition/marginal_sparse.rs#L133)
- `merged_edge.set_fitch_subs(merged_subs)` [partition/marginal_sparse.rs#L195](../../packages/treetime/src/partition/marginal_sparse.rs#L195)

### Swap ref and qry on path edges during reroot

Edges along the old-root-to-new-root path reverse direction during reroot. Ref and qry are swapped on each substitution to match the new parent-to-child direction.

- `apply_reroot()` [partition/marginal_sparse.rs#L133](../../packages/treetime/src/partition/marginal_sparse.rs#L133)
- `edge_data.invert_fitch_subs()` [partition/marginal_sparse.rs#L209](../../packages/treetime/src/partition/marginal_sparse.rs#L209)

### Derive new root sequence from inverted subs

After edge inversion, walk the inverted fitch subs to derive the new root sequence. Each sub's ref state (post-inversion) gives the new-root-ward state at that position.

- `apply_reroot()` [partition/marginal_sparse.rs#L133](../../packages/treetime/src/partition/marginal_sparse.rs#L133)
- `edge_data.fitch_subs()` [partition/marginal_sparse.rs#L235](../../packages/treetime/src/partition/marginal_sparse.rs#L235)

### Reconstruct MAP sequence from parent

Reconstruct a node's full sequence by applying fitch subs and indels to the parent sequence. Used during MAP sequence reconstruction.

- `reconstruct_map_seq()` [partition/marginal_sparse.rs#L90](../../packages/treetime/src/partition/marginal_sparse.rs#L90)
- `edge.fitch_subs()` [partition/marginal_sparse.rs#L98](../../packages/treetime/src/partition/marginal_sparse.rs#L98)

### Count edge mutations for prune decisions

Count total mutations across all partitions for an edge. Used to determine whether a leaf or internal node has enough signal to keep or should be pruned.

- `parse_node_names()` [prune/run.rs#L99](../../packages/treetime/src/commands/prune/run.rs#L99)
- `edge.fitch_subs().len()` [prune/run.rs#L146](../../packages/treetime/src/commands/prune/run.rs#L146)

### Find shared subs between sibling edges

Read fitch subs from both sibling edges and compute their intersection. Shared substitutions will be lifted to the parent edge when siblings are merged.

- `compute_shared_subs_across_partitions()` [topology_cleanup/merge_shared_mutations.rs#L140](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L140)
- `e.fitch_subs()` pair A [topology_cleanup/merge_shared_mutations.rs#L150](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L150)
- `e.fitch_subs()` pair B [topology_cleanup/merge_shared_mutations.rs#L151](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L151)

### Compute private mutations per sibling

Read fitch subs from both sibling edges to compute remaining (non-shared) substitutions. Each child keeps only its private mutations after shared ones move to the parent.

- `merge_sibling_pair()` [topology_cleanup/merge_shared_mutations.rs#L183](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L183)
- `e.fitch_subs()` pair A [topology_cleanup/merge_shared_mutations.rs#L209](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L209)
- `e.fitch_subs()` pair B [topology_cleanup/merge_shared_mutations.rs#L210](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L210)

### Redistribute shared and private subs after sibling merge

After splitting shared vs private mutations, write the shared set onto the new parent edge and the private remainders onto each child edge. Completes the sibling merge operation.

- `merge_sibling_pair()` [topology_cleanup/merge_shared_mutations.rs#L183](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L183)
- `parent_edge.set_fitch_subs()` shared subs [topology_cleanup/merge_shared_mutations.rs#L303](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L303)
- `edge_a.set_fitch_subs()` remaining A [topology_cleanup/merge_shared_mutations.rs#L307](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L307)
- `edge_b.set_fitch_subs()` remaining B [topology_cleanup/merge_shared_mutations.rs#L311](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L311)

### Compose subs when collapsing zero-length edge

When a zero-length internal edge is collapsed, its substitutions are composed with the child edge's substitutions and stored on the surviving child edge.

- `collapse_edge()` [optimize/topology/collapse.rs#L33](../../packages/treetime/src/optimize/topology/collapse.rs#L33)
- `removed_edge_data.chain_fitch_subs(child_edge.fitch_subs())` [optimize/topology/collapse.rs#L60](../../packages/treetime/src/optimize/topology/collapse.rs#L60)
- `child_edge.set_fitch_subs(merged_subs)` [optimize/topology/collapse.rs#L61](../../packages/treetime/src/optimize/topology/collapse.rs#L61)

### Detect mutation-free edges for collapse candidates

Identify mutation-free internal edges as candidates for zero-length collapse. An edge with no fitch subs and no indels across all partitions is considered zero-length.

- `compute_iteration_likelihood()` [optimize/run.rs#L389](../../packages/treetime/src/commands/optimize/run.rs#L389)
- `e.fitch_subs().is_empty()` [optimize/run.rs#L583](../../packages/treetime/src/commands/optimize/run.rs#L583)

### Build likelihood coefficients from variable positions

Collect variable positions from fitch subs (alongside marginal messages) to build the set of positions needing per-edge likelihood coefficients. For each position, look up the fitch sub to determine the parent-to-child state transition for the likelihood computation.

- `get_coefficients()` [partition/optimize_sparse.rs#L45](../../packages/treetime/src/partition/optimize_sparse.rs#L45)
- `edge.fitch_subs().iter().map(Sub::pos)` [partition/optimize_sparse.rs#L58](../../packages/treetime/src/partition/optimize_sparse.rs#L58)
- `edge.fitch_subs().iter().find()` [partition/optimize_sparse.rs#L66](../../packages/treetime/src/partition/optimize_sparse.rs#L66)

### Count mutation pairs for GTR model inference

Count substitution pairs (ref, qry) across all edges to build the mutation count matrix for GTR model inference. Runs before marginal inference, so only fitch subs are available.

- `get_mutation_counts_sparse()` [gtr/infer_gtr/common.rs#L35](../../packages/treetime/src/gtr/infer_gtr/common.rs#L35)
- `e.fitch_subs()` [gtr/infer_gtr/common.rs#L75](../../packages/treetime/src/gtr/infer_gtr/common.rs#L75)
