# Segments with different trees cannot share a run

Segmented genomes such as influenza can have one tree per segment, because reassortment gives the segments different histories. TreeTime takes one `--tree` per run (for example [packages/app-commands/src/commands/shared/tree_input.rs#L58](../../packages/app-commands/src/commands/shared/tree_input.rs#L58)), and every partition of the `--alignment` files lies on that tree ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)). Segments that share one tree run together today; segments with different trees need one run per segment.

## Current state

The partition system attaches several partitions to one graph. [kb/proposals/edge-scoped-inheritance.md](../proposals/edge-scoped-inheritance.md) proposes one graph for all segments: segment trees with named internal nodes merge into one graph, and each edge records which segments it carries. The reader exists as a library function: `segment_trees_read_files()` merges named segment trees by node name into one graph with edge scopes ([kb/algo/graph.md](../algo/graph.md#merging-trees-by-name)). Remaining:

- The command interface: repeated `--tree`, paired with `--alignment` in order, and later the YAML configuration
- Attaching each alignment to the leaves of its segment, and mapping the segment scopes to the partitions of that alignment
- The scope-aware passes ([M-partition-passes-require-tree-topology.md](M-partition-passes-require-tree-topology.md)), without which a merged graph cannot be reconstructed

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

## Proposed implementation plan

The sharing rule is open (Decision required above). The plan stops at that rule.

1. Decide per parameter (clock rate, `Tc` or skyline, substitution model) whether the trees of one run share it. When every parameter is separate, close this issue as "one run per segment" and document the workaround in `docs/en/`, because a run with several trees then adds nothing over a batch of runs
2. When at least one parameter is shared: let `struct TreeArgs` (`packages/app-commands/src/commands/shared/tree_input.rs`) take one tree per segment group, map each `--alignment` file or `--partition` to one tree (config file design in [kb/proposals/config-file-multi-partition.md](../proposals/config-file-multi-partition.md)), and build one graph per group
3. Fit the shared parameters jointly: sum the log-likelihoods of the groups in the objective of each shared parameter (for example the clock regression and the coalescent `Tc` optimization of timetree), and keep the separate parameters per graph. This step is a scientific change and needs a `kb/decisions/` entry with user approval
4. Tests: simulated data from `treetime-simulation` with two trees that share one clock rate; the joint fit recovers the true rate within the interval of the fit, and with separate parameters the run equals two single-tree runs (oracle: the simulation truth, and the single-tree runs)
5. Validate: `just lint-rs`, `just test-rs`, `just generated-check`
