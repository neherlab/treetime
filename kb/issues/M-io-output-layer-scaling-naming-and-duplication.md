# Output layer: quadratic DOT output, ad hoc node names, mixed precision and duplicated code

The code that writes command outputs (`packages/app-output`, the writers in `packages/treetime-io`, and the output steps of `packages/app-commands/src/commands/<command>/run.rs`) has defects in scaling, node naming and number formatting, and structural duplication. Each section states the problem and the required behavior.

## DOT output writes an invisible edge for every pair of leaves

`print_fake_edges()` writes an invisible edge for every ordered pair of roots and of leaves, to keep them in one row ([packages/treetime-io/src/graphviz.rs#L135-L154](../../packages/treetime-io/src/graphviz.rs#L135-L154), called at [#L67](../../packages/treetime-io/src/graphviz.rs#L67) and [#L79](../../packages/treetime-io/src/graphviz.rs#L79)). 2000 leaves give 4 million lines.

Required behavior: stop writing these edges. Delete `print_fake_edges()` and its weight constants.

## Timetree reconstructs sequences for outputs that do not use them

`sequence_outputs_requested()` ([packages/app-commands/src/commands/timetree/run.rs#L473-L480](../../packages/app-commands/src/commands/timetree/run.rs#L473-L480)) enables the final sequence reconstruction for every tree output except graph JSON and DOT. Plain-style Newick and Nexus write no mutations, so a run that asks only for them pays for a reconstruction it discards.

Required behavior: enable the reconstruction only for outputs that read mutations or sequences: Auspice JSON, UShER MAT, Newick and Nexus in an annotated style (`--output-nwk-style=beast` or `nhx`), the reconstructed FASTA, and `--divergence-units=mutations`.

## Unnamed nodes are named differently by each writer

Input parsing names unnamed nodes with `assign_node_names()` ([packages/treetime-graph/src/assign_node_names.rs#L7](../../packages/treetime-graph/src/assign_node_names.rs#L7)), and several pipelines call it again later (timetree, optimize, prune). Nodes created during a run, such as the new root of a `clock` reroot, can stay unnamed, and each writer then invents its own name: an empty string in UShER MAT ([packages/app-output/src/usher_mat.rs#L90](../../packages/app-output/src/usher_mat.rs#L90)), `node_<key>` in Auspice JSON, augur node data, the mugration CSVs and the amino-acid FASTA (`node_name_or_key()`, [packages/treetime-graph/src/assign_node_names.rs#L55-L57](../../packages/treetime-graph/src/assign_node_names.rs#L55-L57)), no label in Newick, and `NODE_0000001` in the reconstructed FASTA. The same node gets different names in the files of one run.

Required behavior: the project convention is `NODE_0000001` (7 digits, numbered in preorder, skipping names already in use), as v0 names nodes in `_prepare_nodes()` ([packages/legacy/treetime/treetime/treeanc.py#L460-L479](../../packages/legacy/treetime/treetime/treeanc.py#L460-L479)). Naming happens outside the formats: during input parsing for nodes that exist in the input, and at node creation for nodes created during a run (new root, polytomy resolution, merged nodes). Every node then has a name, the name map holds `String` instead of `Option<String>`, the repeated `assign_node_names()` calls go away, and writers never invent names.

## Mugration probabilities are written with three precisions

One per-node state profile is written three ways: the confidence CSV writes every state with 6 decimals ([packages/app-output/src/trait_tables.rs#L56](../../packages/app-output/src/trait_tables.rs#L56)), augur node data writes states above 0.001 at full precision, and Auspice JSON writes states above 0.001 rounded to 3 significant digits, and the entropy too ([packages/app-output/src/auspice.rs#L315-L324](../../packages/app-output/src/auspice.rs#L315-L324); the 0.001 filter is `build_confidence_map()`, [packages/app-output/src/trait_profile.rs#L6-L11](../../packages/app-output/src/trait_profile.rs#L6-L11)).

Required behavior: every output writes probabilities and the entropy at full precision. The state filters stay as they are.

## Graph traversal allocates on every visit

`GraphNodeForward::new()` and `GraphNodeBackward::new()` collect the parent and child keys of each visited node into new `Vec`s ([packages/treetime-graph/src/graph_traverse.rs#L131-L137](../../packages/treetime-graph/src/graph_traverse.rs#L131-L137), [#L160-L166](../../packages/treetime-graph/src/graph_traverse.rs#L160-L166)), two allocations per node per traversal.

Required behavior: no per-visit allocation; callers that need the neighbor keys read them from the graph.

## The output vocabulary is spelled out in six places

Every output is listed separately in `OutputSelection` and in `CommandKind::non_tree_outputs()` and `default_outputs()` ([packages/app-output/src/output_plan.rs](../../packages/app-output/src/output_plan.rs)), in the per-command selection enums made by `per_command_output_selection!` ([packages/app-commands/src/commands/shared/output_args.rs](../../packages/app-commands/src/commands/shared/output_args.rs)), in the per-file flag lists of `impl_resolve_outputs!` ([packages/app-commands/src/commands/shared/resolve_outputs.rs](../../packages/app-commands/src/commands/shared/resolve_outputs.rs)), in the per-file flag fields of each command's `args.rs`, and in `result_outputs()` that the UI reader uses ([packages/app-commands/src/results/outputs.rs](../../packages/app-commands/src/results/outputs.rs)). `flag_name()` restates the clap flag names. Adding one output means editing every list, and the lists can drift apart: the clock chart outputs (`clock-chart-svg`, `clock-chart-png`) needed an entry in each list except `result_outputs()`.

Required behavior: one table with one row per output (name, file extension, CLI flag, description, commands, default status) that every list derives from.


## Mugration duplicates ancestral reconstruction

Mugration is a strict subset of ancestral reconstruction: the reconstruction of one discrete character on the tree. v0 implements it exactly that way ([packages/legacy/treetime/treetime/wrappers.py#L699-L800](../../packages/legacy/treetime/treetime/wrappers.py#L699-L800)): it maps each state to one letter of a custom alphabet and missing data to an ambiguous letter, builds a one-column alignment of the leaves, and runs `TreeAnc.infer_ancestral_sequences(method='ml', marginal=True, infer_gtr=True, reconstruct_tip_states=True)` with a custom GTR (all-ones `W`, equilibrium frequencies from the optional weights), then `optimize_gtr_rate()`.

v1 reimplements it beside ancestral: its own pipeline ([packages/treetime/src/mugration/pipeline.rs](../../packages/treetime/src/mugration/pipeline.rs)), its own partition type ([packages/treetime/src/partition/marginal/discrete/](../../packages/treetime/src/partition/marginal/discrete/)), and its own output code for the state confidence map, the entropy and the transition labels (`trait_profile.rs`, `trait_tables.rs`, `augur_node_data_traits.rs` and the trait part of `auspice.rs` in `packages/app-output/src/`). This adds ad hoc code, a parallel structure and duplicated computations that drift from ancestral.

Required behavior: reimplement mugration on top of ancestral reconstruction. Mugration-specific code reduces to building the state alphabet and the one-site alignment from the metadata and mapping reconstructed letters back to states; reconstruction, GTR inference, confidence and entropy come from the ancestral code. The ancestral `Alphabet` ([packages/treetime/src/alphabet/alphabet.rs](../../packages/treetime/src/alphabet/alphabet.rs)) is built only from named alphabets (`AlphabetName`: nucleotide, amino acid) and needs one built from the trait states. Delete the discrete partition and the duplicated computations afterwards.

## KB entries describe removed code

- `kb/decisions/multi-format-tree-io.md` describes format traits, `ConverterGraph` and a `convert` command that no longer exist
- `kb/issues/H-graph-lock-topology-is-exposed-to-consumers.md` describes `Arc<RwLock<Node>>`; the graph has no locks now
- `kb/issues/M-command-output-ownership-is-scattered.md` and `kb/issues/M-output-module-mixes-topology-ordering.md` cite paths under `packages/treetime/src/commands/`, which no longer exists

Required behavior: update each entry to the current code, or delete it when the problem it records is gone.

## Intended behavior

Each output file is written as soon as its data is ready. A run that fails late keeps the files it finished, so users can salvage part of the results.

## Related issues

- [N-clock-unnamed-root-after-reroot.md](N-clock-unnamed-root-after-reroot.md)
- [N-io-confidence-tsv-unnamed-root-has-empty-name.md](N-io-confidence-tsv-unnamed-root-has-empty-name.md)
- [N-mugration-confidence-rows-copied-for-output.md](N-mugration-confidence-rows-copied-for-output.md)
- [N-kb-multi-format-tree-io-decision-describes-removed-converters.md](N-kb-multi-format-tree-io-decision-describes-removed-converters.md)
