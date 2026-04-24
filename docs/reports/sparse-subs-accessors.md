# Sparse substitution accessors

Production call sites for `marginal_subs` and `fitch_subs` on `SparseEdgePartition`, excluding tests.

Definitions: [packages/treetime/src/representation/payload/sparse.rs](../../packages/treetime/src/representation/payload/sparse.rs)

## Marginal subs writes (5)

### Store computed marginal subs after forward pass

After the child profile is finalized during the forward (root-to-tip) marginal pass, compute the marginal substitutions for each parent edge by comparing parent and child marginal profiles, then store them on the edge.

- `process_node_forward()` [representation/partition/marginal_passes.rs#L250](../../packages/treetime/src/representation/partition/marginal_passes.rs#L250)
- `compute_marginal_subs_for_edge()` [representation/partition/marginal_passes.rs#L339](../../packages/treetime/src/representation/partition/marginal_passes.rs#L339)
- `edge.set_marginal_subs(marginal_subs)` [representation/partition/marginal_passes.rs#L341](../../packages/treetime/src/representation/partition/marginal_passes.rs#L341)

### Invalidate marginal subs on edge split during reroot

The old edge moves to the child side of a new internal node. Marginal subs are cleared because they were computed under the old root and are stale.

- `apply_reroot()` [representation/partition/marginal_sparse.rs#L133](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L133)
- `old_edge_data.clear_marginal_subs()` [representation/partition/marginal_sparse.rs#L146](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L146)

### Invalidate marginal subs after edge inversion during reroot

After inverting fitch subs on edges along the old-root-to-new-root path, marginal subs are cleared because they depend on directionality that just changed.

- `apply_reroot()` [representation/partition/marginal_sparse.rs#L133](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L133)
- `edge_data.clear_marginal_subs()` [representation/partition/marginal_sparse.rs#L212](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L212)

## Marginal subs reads (3)

### Count marginal subs for initial branch length guess

Count discrete substitutions per edge across all partitions to estimate initial branch lengths before the optimization loop. Called from `optimize` and `timetree` commands.

- `initial_guess_mixed()` [commands/optimize/optimize_unified.rs#L729](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L729)
- `partition.edge_subs()` [commands/optimize/optimize_unified.rs#L774](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L774)
- `edge.marginal_subs()` [representation/partition/marginal_sparse.rs#L269](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L269)

### Collect marginal subs for branch mutation annotation

Collect all marginal subs across partitions for each edge and write them as comma-separated mutation strings onto the graph nodes. Used for output in `ancestral`, `optimize`, and `timetree` commands.

- `annotate_branch_mutations()` [representation/payload/ancestral.rs#L152](../../packages/treetime/src/representation/payload/ancestral.rs#L152)
- `partition.edge_subs()` [representation/payload/ancestral.rs#L170](../../packages/treetime/src/representation/payload/ancestral.rs#L170)
- `edge.marginal_subs()` [representation/partition/marginal_sparse.rs#L269](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L269)

## Fitch subs call sites (22)

### Append resolved subs during Fitch forward pass

During the forward (tip-to-root) Fitch parsimony pass, substitutions are resolved by comparing child and parent Fitch sets at each variable position. Resolved subs are appended to the edge.

- `run_fitch_forward()` [commands/ancestral/fitch.rs#L292](../../packages/treetime/src/commands/ancestral/fitch.rs#L292)
- `edge.extend_fitch_subs(subs)` [commands/ancestral/fitch.rs#L453](../../packages/treetime/src/commands/ancestral/fitch.rs#L453)

### Apply fitch subs to reconstruct child sequence for output

After Fitch reconstruction, apply fitch subs to the parent sequence to produce the child node's reconstructed sequence for output.

- `run_fitch_reconstruction()` [commands/ancestral/fitch.rs#L561](../../packages/treetime/src/commands/ancestral/fitch.rs#L561)
- `edge_part.fitch_subs()` [commands/ancestral/fitch.rs#L580](../../packages/treetime/src/commands/ancestral/fitch.rs#L580)

### Seed candidate positions for marginal sub computation

Fitch subs provide candidate positions (union with parent/child variable sites) where marginal substitutions may differ from Fitch. The marginal computation refines these using posterior probabilities.

- `compute_marginal_subs_for_edge()` [representation/partition/marginal_passes.rs#L81](../../packages/treetime/src/representation/partition/marginal_passes.rs#L81)
- `.fitch_subs()` [representation/partition/marginal_passes.rs#L100](../../packages/treetime/src/representation/partition/marginal_passes.rs#L100)

### Identify variable positions in backward marginal pass

During the backward (tip-to-root) marginal pass, fitch subs identify variable positions on each child edge. These positions seed the parent's variable-site map for marginal profile computation.

- `process_node_backward()` [representation/partition/marginal_passes.rs#L128](../../packages/treetime/src/representation/partition/marginal_passes.rs#L128)
- `edge_data.fitch_subs()` [representation/partition/marginal_passes.rs#L174](../../packages/treetime/src/representation/partition/marginal_passes.rs#L174)

### Populate state maps in forward marginal pass

During the forward (root-to-tip) marginal pass, fitch subs on parent edges identify positions where the parent and child states differ. These positions populate the parent/child state maps for marginal profile refinement.

- `process_node_forward()` [representation/partition/marginal_passes.rs#L250](../../packages/treetime/src/representation/partition/marginal_passes.rs#L250)
- `edge_data.fitch_subs()` [representation/partition/marginal_passes.rs#L272](../../packages/treetime/src/representation/partition/marginal_passes.rs#L272)

### Provide state pairs for child message precalculation

During message precalculation in the forward pass, fitch subs on each child edge provide parent/child state pairs. These are combined with sibling messages to compute the outgoing message to each child.

- `process_node_forward()` [representation/partition/marginal_passes.rs#L250](../../packages/treetime/src/representation/partition/marginal_passes.rs#L250)
- `child_edge_data.fitch_subs()` [representation/partition/marginal_passes.rs#L364](../../packages/treetime/src/representation/partition/marginal_passes.rs#L364)

### Update root sequence from parent edge during reroot merge

During edge merge in reroot, fitch subs on the parent edge update the root sequence to the old root's state (applying ref states at each sub position).

- `apply_reroot()` [representation/partition/marginal_sparse.rs#L133](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L133)
- `edge.fitch_subs()` [representation/partition/marginal_sparse.rs#L171](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L171)

### Compose parent and child subs during reroot edge merge

When reroot removes an internal node, the parent and child edges are merged into one. The two substitution lists are composed into a single consistent list for the merged edge.

- `apply_reroot()` [representation/partition/marginal_sparse.rs#L133](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L133)
- `parent_edge.chain_fitch_subs(child_edge.fitch_subs())` [representation/partition/marginal_sparse.rs#L189](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L189)

### Store composed subs on merged edge after reroot

Store the composed substitution list on the newly created merged edge after edge merge.

- `apply_reroot()` [representation/partition/marginal_sparse.rs#L133](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L133)
- `merged_edge.set_fitch_subs(merged_subs)` [representation/partition/marginal_sparse.rs#L195](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L195)

### Swap ref and qry on path edges during reroot

Edges along the old-root-to-new-root path reverse direction during reroot. Ref and qry are swapped on each substitution to match the new parent-to-child direction.

- `apply_reroot()` [representation/partition/marginal_sparse.rs#L133](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L133)
- `edge_data.invert_fitch_subs()` [representation/partition/marginal_sparse.rs#L209](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L209)

### Derive new root sequence from inverted subs

After edge inversion, walk the inverted fitch subs to derive the new root sequence. Each sub's ref state (post-inversion) gives the new-root-ward state at that position.

- `apply_reroot()` [representation/partition/marginal_sparse.rs#L133](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L133)
- `edge_data.fitch_subs()` [representation/partition/marginal_sparse.rs#L235](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L235)

### Reconstruct MAP sequence from parent

Reconstruct a node's full sequence by applying fitch subs and indels to the parent sequence. Used during MAP sequence reconstruction.

- `reconstruct_map_seq()` [representation/partition/marginal_sparse.rs#L90](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L90)
- `edge.fitch_subs()` [representation/partition/marginal_sparse.rs#L98](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L98)

### Count edge mutations for prune decisions

Count total mutations across all partitions for an edge. Used to determine whether a leaf or internal node has enough signal to keep or should be pruned.

- `parse_node_names()` [commands/prune/run.rs#L120](../../packages/treetime/src/commands/prune/run.rs#L120)
- `edge.fitch_subs().len()` [commands/prune/run.rs#L167](../../packages/treetime/src/commands/prune/run.rs#L167)

### Find shared subs between sibling edges

Read fitch subs from both sibling edges and compute their intersection. Shared substitutions will be lifted to the parent edge when siblings are merged.

- `compute_shared_subs_across_partitions()` [commands/prune/run.rs#L406](../../packages/treetime/src/commands/prune/run.rs#L406)
- `e.fitch_subs()` pair A [commands/prune/run.rs#L416](../../packages/treetime/src/commands/prune/run.rs#L416)
- `e.fitch_subs()` pair B [commands/prune/run.rs#L417](../../packages/treetime/src/commands/prune/run.rs#L417)

### Compute private mutations per sibling

Read fitch subs from both sibling edges to compute remaining (non-shared) substitutions. Each child keeps only its private mutations after shared ones move to the parent.

- `merge_sibling_pair()` [commands/prune/run.rs#L449](../../packages/treetime/src/commands/prune/run.rs#L449)
- `e.fitch_subs()` pair A [commands/prune/run.rs#L475](../../packages/treetime/src/commands/prune/run.rs#L475)
- `e.fitch_subs()` pair B [commands/prune/run.rs#L476](../../packages/treetime/src/commands/prune/run.rs#L476)

### Redistribute shared and private subs after sibling merge

After splitting shared vs private mutations, write the shared set onto the new parent edge and the private remainders onto each child edge. Completes the sibling merge operation.

- `merge_sibling_pair()` [commands/prune/run.rs#L449](../../packages/treetime/src/commands/prune/run.rs#L449)
- `parent_edge.set_fitch_subs()` shared subs [commands/prune/run.rs#L569](../../packages/treetime/src/commands/prune/run.rs#L569)
- `edge_a.set_fitch_subs()` remaining A [commands/prune/run.rs#L573](../../packages/treetime/src/commands/prune/run.rs#L573)
- `edge_b.set_fitch_subs()` remaining B [commands/prune/run.rs#L577](../../packages/treetime/src/commands/prune/run.rs#L577)

### Compose subs when collapsing zero-length edge

When a zero-length internal edge is collapsed, its substitutions are composed with the child edge's substitutions and stored on the surviving child edge.

- `collapse_edge()` [representation/algo/topology_cleanup/collapse.rs#L33](../../packages/treetime/src/representation/algo/topology_cleanup/collapse.rs#L33)
- `removed_edge_data.chain_fitch_subs(child_edge.fitch_subs())` [representation/algo/topology_cleanup/collapse.rs#L60](../../packages/treetime/src/representation/algo/topology_cleanup/collapse.rs#L60)
- `child_edge.set_fitch_subs(merged_subs)` [representation/algo/topology_cleanup/collapse.rs#L61](../../packages/treetime/src/representation/algo/topology_cleanup/collapse.rs#L61)

### Detect mutation-free edges for collapse candidates

Identify mutation-free internal edges as candidates for zero-length collapse. An edge with no fitch subs and no indels across all partitions is considered zero-length.

- `compute_iteration_likelihood()` [commands/optimize/run.rs#L389](../../packages/treetime/src/commands/optimize/run.rs#L389)
- `e.fitch_subs().is_empty()` [commands/optimize/run.rs#L583](../../packages/treetime/src/commands/optimize/run.rs#L583)

### Build likelihood coefficients from variable positions

Collect variable positions from fitch subs (alongside marginal messages) to build the set of positions needing per-edge likelihood coefficients. For each position, look up the fitch sub to determine the parent-to-child state transition for the likelihood computation.

- `get_coefficients()` [commands/optimize/optimize_sparse.rs#L45](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L45)
- `edge.fitch_subs().iter().map(Sub::pos)` [commands/optimize/optimize_sparse.rs#L58](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L58)
- `edge.fitch_subs().iter().find()` [commands/optimize/optimize_sparse.rs#L66](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L66)

### Count mutation pairs for GTR model inference

Count substitution pairs (ref, qry) across all edges to build the mutation count matrix for GTR model inference. Runs before marginal inference, so only fitch subs are available.

- `get_mutation_counts_sparse()` [gtr/infer_gtr/sparse.rs#L35](../../packages/treetime/src/gtr/infer_gtr/sparse.rs#L35)
- `e.fitch_subs()` [gtr/infer_gtr/sparse.rs#L75](../../packages/treetime/src/gtr/infer_gtr/sparse.rs#L75)
