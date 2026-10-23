# Segments with different trees cannot share a run

Segmented genomes such as influenza can have one tree per segment, because reassortment gives the segments different histories. TreeTime takes one `--tree` per run (for example [packages/app-commands/src/commands/timetree/args.rs#L129](../../packages/app-commands/src/commands/timetree/args.rs#L129)), and every partition of the `--alignment` files lies on that tree ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)). Segments that share one tree run together today; segments with different trees need one run per segment.

## Current state

The partition system attaches several partitions to one graph. A run with one tree per segment would need one graph per segment group and a rule for the parameters that the groups share.

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
