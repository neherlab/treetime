# Runs over several segment trees cover three commands and lack a rule for shared parameters

Segmented genomes such as influenza can have one tree per segment, because reassortment gives the segments different histories. TreeTime takes one `--tree` per run (for example [packages/app-commands/src/commands/shared/tree_input.rs#L58](../../packages/app-commands/src/commands/shared/tree_input.rs#L58)), and every partition of the `--alignment` files lies on that tree ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)). `--segment-trees` runs `ancestral`, `optimize` and `timetree` on all segment trees at once (current state below). The other commands, `--translations`, and the rule for the parameters that the segments share are open.

## Current state

The partition system attaches several partitions to one graph. [kb/proposals/edge-scoped-inheritance.md](../proposals/edge-scoped-inheritance.md) proposes one graph for all segments: segment trees with named internal nodes merge into one graph, and each edge records which segments it carries. `treetime ancestral`, `treetime optimize` and `treetime timetree` take `--segment-trees <tree>... --segment-alignments <fasta>...`: one tree and one alignment per segment, merged by node name into one graph with edge scopes ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name)). `ancestral` reconstructs each segment on that graph by Fitch parsimony or marginal likelihood ([kb/algo/ancestral.md](../algo/ancestral.md)); `optimize` fits one length per shared edge ([kb/algo/optimization.md](../algo/optimization.md#branch-lengths-on-graphs-of-segment-trees)). `timetree` dates the graph at one clock rate for all segments, given with `--clock-rate` ([kb/issues/M-timetree-segment-graphs-fixed-clock-only.md](M-timetree-segment-graphs-fixed-clock-only.md)); the segments have their own relative rates, as in optimize. The flags replace `--tree` and `--alignment` and pair trees with alignments in order; they are separate from `--tree` so that the configuration of the other commands is unchanged.

Outputs: the per-alignment outputs (reconstructed sequences, Nexus, Auspice JSON, annotated Newick, UShER MAT, augur node data of `ancestral`) hold the tree of each segment, cut from the graph by its scopes, with nodes that the segment passes through as nodes with one child. Graph JSON and DOT hold the whole graph. Outputs that hold one tree for the whole run (plain Newick without `{alignment}` in its path, augur node data of `optimize`) are skipped with a warning. Remaining:

- `--translations` with `--segment-trees` (an error now)
- Graph JSON holds topology only: node names and edge scopes are not in it
- The other commands, and one interface with `--tree`, `--alignment`, and the YAML configuration

Toy input for this reader is in [data/toy/reassortment/](../../data/toy/reassortment/). Converting TreeKnit output (local internal names, pairwise MCCs) into named trees needs a prototype script: [M-io-treeknit-output-reaches-segment-graphs-only-through-a-python-prototype.md](M-io-treeknit-output-reaches-segment-graphs-only-through-a-python-prototype.md).

> [!IMPORTANT]
> **Decision required.** A segment graph shares its edges and node dates by design ([kb/decisions/segment-graph-replaces-v0-arg-masks.md](../decisions/segment-graph-replaces-v0-arg-masks.md)), `timetree` takes one clock rate for all segments, and each segment has its own substitution model and relative rate. No rule defines the other parameters of timetree: an estimated clock rate and the coalescent model (the merger rate `Tc` or the skyline). Options per parameter:
>
> - Shared: one value fitted jointly over all segments of the run
> - Separate: one value per segment

The public alignment-input design discussion proposes separate nucleotide and amino-acid alignment partitions [[issue](https://github.com/neherlab/treetime/issues/306)]. A maintainer response identifies sequence alphabets, discrete traits, sparse/dense storage, and configuration-based partition specifications as parts of the same partition abstraction [[comment](https://github.com/neherlab/treetime/issues/306#issuecomment-2565383822)].

## v0 comparison

v0 accepts one tree and one alignment file. Users run each segment on its own tree, or concatenate segments that share a tree into one alignment.

## Workaround

For the other commands and for `--translations`: run TreeTime once per segment with the tree of that segment.

## Related

### Known issues

- [N-io-sequence-name-matching-unreliable.md](N-io-sequence-name-matching-unreliable.md) -- name matching affects multi-segment attachment
- [N-io-large-dataset-memory-constraint.md](N-io-large-dataset-memory-constraint.md) -- memory constraints compound with multiple segments

### Proposals

- [kb/proposals/config-file-multi-partition.md](../proposals/config-file-multi-partition.md) -- configuration file format that would subsume per-segment CLI flags
- [kb/proposals/unified-input-format-support.md](../proposals/unified-input-format-support.md) -- alternative input via formats with embedded sequences

## Proposed implementation plan

The sharing rule is open (Decision required above). The plan stops at that rule.

1. Decide per parameter (estimated clock rate, `Tc` or skyline) whether the segments of one run share it
2. Fit each parameter on the segment graph, in the clock rate estimation and the coalescent prior of `timetree` on segment graphs ([M-timetree-segment-graphs-fixed-clock-only.md](M-timetree-segment-graphs-fixed-clock-only.md)): a shared parameter sums the log-likelihoods of the segments in its objective, a separate one is fitted per segment. This step is a scientific change and needs a `kb/decisions/` entry with user approval
3. Tests: simulated segment graphs from `treetime-simulation` whose segments share one clock rate; the joint fit recovers the true rate within the interval of the fit (oracle: the simulation truth)
4. Validate: `just lint-rs`, `just test-rs`, `just generated-check`
