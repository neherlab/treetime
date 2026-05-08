# Standalone tree pruning command

v1 exposes tree pruning as a standalone CLI command `treetime prune` with three independent pruning criteria: short branches, empty branches (zero mutations), and explicit node names. v0 has only an internal `prune_short_branches()` method with a hardcoded probabilistic threshold, not accessible from the CLI.

## What v0 does

v0 prunes short branches inside `TreeAnc.prune_short_branches()` ([packages/legacy/treetime/treetime/treeanc.py#L1475-L1496](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1496)), called automatically during ancestral reconstruction when `prune_short=True` (the default). The method operates on internal nodes only and applies two criteria together:

1. Branch length below `0.1 * one_mutation`, where `one_mutation = 1.0 / sequence_length` (relative threshold scaled by alignment size).
2. Probability of identical parent-child sequences exceeds 0.1, computed via `gtr.prob_t(edge.branch_length)` (GTR transition probability at zero divergence).

Both conditions must hold for pruning to occur. The threshold is not user-configurable.

v0 does not support pruning by mutation count or by explicit node name. There is a `--prune-outliers` flag in the timetree argument parser, but that removes temporal outliers (sequences whose dates are inconsistent with the molecular clock), not topologically short or empty branches.

## What v1 does

v1 adds a dedicated `treetime prune` subcommand ([packages/treetime/src/commands/prune/run.rs#L28-L81](../../packages/treetime/src/commands/prune/run.rs#L28-L81)) that operates on an input tree (and optionally an alignment) and writes a pruned tree.

### CLI arguments

| Argument                  | Default    | Purpose                                                     |
| ------------------------- | ---------- | ----------------------------------------------------------- |
| `--tree`                  | (required) | Input tree in Newick/Nexus/Phylip                           |
| `--aln`                   | (none)     | FASTA alignment (required if `--prune-empty`)               |
| `--alphabet`              | auto       | Nucleotide or amino acid alphabet                           |
| `--outdir`                | (required) | Output directory                                            |
| `--prune-short`           | (none)     | Threshold for short branch pruning (substitutions per site) |
| `--prune-empty`           | false      | Prune branches with zero mapped mutations                   |
| `--prune-nodes-list`      | (none)     | Comma-separated node names to remove                        |
| `--prune-nodes-list-file` | (none)     | File with node names to remove (one per line)               |

### Algorithm

The command runs two passes over the tree:

Pass 1: Internal node pruning ([packages/treetime/src/commands/prune/run.rs#L154-L197](../../packages/treetime/src/commands/prune/run.rs#L154-L197)). For each edge targeting an internal node, any of three independent criteria triggers pruning:

- Short branch: branch length < user-specified threshold (strict less-than).
- Empty branch: mutation data is present for the edge AND total mutation count is zero across all partitions. Edges with no partition data (`None`) are preserved: unknown is not the same as zero.
- Named node: target node's name appears in the user-provided set.

When an internal node is pruned, its children are reconnected to its parent via `Graph.collapse_edge()` ([packages/treetime-graph/src/graph_ops.rs#L146-L206](../../packages/treetime-graph/src/graph_ops.rs#L146-L206)). Branch lengths are summed when both edges have values. Substitution lists from the removed edge are merged into each reconnected edge using sorted union with deduplication.

Pass 2: Leaf pruning ([packages/treetime/src/commands/prune/run.rs#L199-L228](../../packages/treetime/src/commands/prune/run.rs#L199-L228)). Only the named-node criterion applies to leaves. Short-branch and empty-branch pruning do not affect leaves. When a leaf is removed, the algorithm recursively collapses parent nodes that become childless or unary (have at most one child and are not root).

Output: `pruned_tree.nwk` and `pruned_tree.nexus` in the output directory.

## Scientific background

### Short branches in phylogenetic trees

Branch lengths in phylogenetic trees represent evolutionary distance in substitutions per site. Very short branches indicate negligible observed divergence between ancestor and descendant. These arise from:

- Insufficient phylogenetic signal. When an alignment is short relative to the divergence time, many true branching events produce zero or near-zero observed substitutions. The branching order among these nodes is effectively unresolved.
- Rapid radiation. When lineages diverge in quick succession (adaptive radiation, founder events, epidemic expansion), the intervening internal branches are genuinely short, and the tree topology at these nodes is a soft polytomy.
- Overresolved input trees. Parsimony and distance methods can produce fully resolved binary trees even when data support only partial resolution. The spurious internal branches have near-zero lengths.

Collapsing short branches converts the resolved but unsupported binary splits into multifurcations (polytomies) that honestly represent the phylogenetic uncertainty. This is standard practice before downstream analysis: bootstrap support thresholds, Bayesian posterior probability cutoffs, and branch length thresholds all serve the same purpose of removing resolution that the data cannot support.

### Zero-mutation branches

A branch with zero mapped substitutions after Fitch parsimony compression carries no phylogenetic information beyond topology. All sequences in the subtree below such a branch are identical at variable positions. These branches arise when:

- Multiple identical sequences exist in the alignment (common in closely related viral isolates during an outbreak).
- The internal branching order among identical sequences was resolved arbitrarily by the tree-building algorithm.

Pruning zero-mutation branches collapses these arbitrary resolutions. The resulting polytomy accurately represents the data: the sequences are indistinguishable at the resolution of the alignment, and their branching order is unknown.

### Taxon pruning and tree surgery

Removing specific taxa (leaves) from a phylogenetic tree is a standard preprocessing step in phylogenetics. Common use cases include removing outlier sequences (contaminated, misidentified, or low-quality), focusing analysis on a geographic or temporal subset, and pruning duplicate or redundant sequences. After leaf removal, the resulting degree-2 internal nodes (a parent with a single child) are collapsed by summing the two branch lengths, preserving the total evolutionary distance.

## Differences from v0's internal pruning

| Aspect                  | v0 (`prune_short_branches`)                           | v1 (`prune` command)                                             |
| ----------------------- | ----------------------------------------------------- | ---------------------------------------------------------------- |
| Exposure                | Internal method, called automatically                 | Standalone CLI subcommand                                        |
| Short branch threshold  | Hardcoded: `0.1 / sequence_length`                    | User-specified via `--prune-short`                               |
| Probabilistic criterion | Requires P(identical sequences) > 0.1 via GTR         | No probability check                                             |
| Empty branch criterion  | Not available                                         | `--prune-empty` with Fitch mutation counts                       |
| Named node removal      | Not available                                         | `--prune-nodes-list`, `--prune-nodes-list-file`                  |
| Leaf pruning            | Not available                                         | By explicit name                                                 |
| Branch length handling  | Children reparented without branch length summation   | Branch lengths summed when both edges have values                |
| Mutation handling       | No mutation tracking                                  | Substitution lists merged via sorted union                       |
| Recursive cleanup       | Not needed (only internal, never creates unary nodes) | Recursive collapse of childless/unary parents after leaf removal |

## Practical impact

Users can run `treetime prune --tree=tree.nwk --aln=aln.fasta --prune-short=1e-6 --prune-empty --outdir=out/` to simplify a tree before downstream analysis. The command is composable with other v1 commands: prune first, then run ancestral reconstruction or timetree inference on the simplified tree.

The three pruning criteria are independent and can be combined freely. Using `--prune-empty` without `--aln` produces an error. Using `--prune-short` or `--prune-nodes-list` without `--aln` works on topology and branch lengths alone.
