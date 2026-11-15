# Segments with different trees cannot share a run

Segmented genomes such as influenza can have one tree per segment, because reassortment gives the segments different histories. TreeTime takes one `--tree` per run (for example [packages/app-commands/src/commands/timetree/args.rs#L129](../../packages/app-commands/src/commands/timetree/args.rs#L129)), and every partition of the `--alignment` files lies on that tree ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)). Segments that share one tree run together today; segments with different trees need one run per segment.

## Current state

The partition system attaches several partitions to one graph. [kb/proposals/edge-scoped-inheritance.md](../proposals/edge-scoped-inheritance.md) proposes one graph for all segments: segment trees with named internal nodes merge into one graph, and each edge records which segments it carries. `treetime ancestral --segment-trees <tree>... --segment-alignments <fasta>...` reads one tree and one alignment per segment, merges the trees by node name into one graph with edge scopes ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name)), and reconstructs each segment on that graph by Fitch parsimony or marginal likelihood, dense or sparse ([kb/algo/ancestral.md](../algo/ancestral.md)). The flags replace `--tree` and `--alignment` and pair trees with alignments in order; they are separate from `--tree` so that the configuration of the other commands is unchanged. Such a run writes the reconstructed sequences per segment, the GTR model, graph JSON, and DOT; it skips the outputs that hold one tree (Newick, Nexus, Auspice JSON, UShER MAT, augur node data) with a warning. Remaining:

- `--translations` with `--segment-trees` (an error now)
- Tree outputs per segment: the tree that each segment's scopes trace, written in the tree formats
- Graph JSON holds topology only: node names and edge scopes are not in it
- The other commands, and one interface with `--tree`, `--alignment`, and the YAML configuration

Toy input for this reader is in [data/toy/reassortment/](../../data/toy/reassortment/). Converting TreeKnit output (local internal names, pairwise MCCs) into named trees is a further step described in the proposal.

> [!IMPORTANT]
> **Decision required.** No rule defines which parameters the trees of one run share. Examples in timetree are the clock rate and the coalescent model (the merger rate `Tc` or the skyline). Options per parameter:
>
> - Shared: one value fitted jointly over all trees of the run
> - Separate: one value per tree, as in one run per segment today
>
> The choice also decides whether a run with several trees is more than a batch of independent runs.

The public alignment-input design discussion proposes separate nucleotide and amino-acid alignment partitions [[issue](https://github.com/neherlab/treetime/issues/306)]. A maintainer response identifies sequence alphabets, discrete traits, sparse/dense storage, and configuration-based partition specifications as parts of the same partition abstraction [[comment](https://github.com/neherlab/treetime/issues/306#issuecomment-2565383822)].

## v0 comparison

v0 accepts one tree and one alignment file. Users run each segment on its own tree, or concatenate segments that share a tree into one alignment.

## Workaround

Run TreeTime once per segment with the tree of that segment.

## Related

### Known issues

- [N-io-sequence-name-matching-unreliable.md](N-io-sequence-name-matching-unreliable.md) -- name matching affects multi-segment attachment
- [N-io-large-dataset-memory-constraint.md](N-io-large-dataset-memory-constraint.md) -- memory constraints compound with multiple segments

### Proposals

- [kb/proposals/config-file-multi-partition.md](../proposals/config-file-multi-partition.md) -- configuration file format that would subsume per-segment CLI flags
- [kb/proposals/unified-input-format-support.md](../proposals/unified-input-format-support.md) -- alternative input via formats with embedded sequences
