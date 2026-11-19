# Influenza A/H3N2: four segment trees

Seasonal influenza A/H3N2 viruses sampled from 2004 to 2026, with one tree and one alignment for each of the segments HA, NA, PB1 and PB2. The segments of influenza reassort: two viruses that infect one cell can exchange segments, so each segment has its own tree. TreeTime merges the four trees into one graph by node name. An ancestor that several segments share has the same name in their trees and is one node of the graph, which therefore gets one date and one set of mutations from all of these segments.

## Datasets

| Folder           | Viruses | Internal nodes (HA, NA, PB1, PB2) | Nodes of the merged graph | Samples before 2021 |
| ---------------- | ------: | --------------------------------- | ------------------------: | ------------------: |
| [`100/`](100/)   |     100 | 65, 51, 54, 64                    |                       266 |                   6 |
| [`500/`](500/)   |     500 | 270, 249, 262, 267                |                      1026 |                  18 |
| [`1564/`](1564/) |    1564 | 778, 747, 767, 783                |                      2839 |                  49 |

Each folder holds:

- `tree_ha.nwk`, `tree_na.nwk`, `tree_pb1.nwk`, `tree_pb2.nwk`: rooted trees in substitutions per site, with a name on every node. A name that occurs in several trees is one ancestor of these segments
- `aln_ha.fasta.xz`, `aln_na.fasta.xz`, `aln_pb1.fasta.xz`, `aln_pb2.fasta.xz`: nucleotide alignments of the leaves (1737, 1436, 2341 and 2341 columns)
- `metadata.tsv`: one row per virus, keyed by `strain`, with the sampling `date` (`2004-XX-XX` when only the year is known), `region`, `country`, the Nextclade subclades of HA and NA, the INSDC accession of each segment, and the submitters

The viruses come from all ten regions of the source; Europe gives about two thirds of them. The oldest samples date the deep nodes of the trees: every set contains A/Wellington/1/2004, A/California/7/2004 and A/Nepal/921/2006, because the subsets take the oldest leaf of each region first.

## Examples

| File                                                                                                                                                                                             | Command                                                               | Demonstrates                                                                                                                   |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| [`100/ancestral.yaml`](100/ancestral.yaml), [`1564/ancestral.yaml`](1564/ancestral.yaml)                                                                                                         | `ancestral`                                                           | Marginal ancestral sequences of each segment on the graph                                                                      |
| [`100/ancestral-parsimony.yaml`](100/ancestral-parsimony.yaml), [`1564/ancestral-parsimony.yaml`](1564/ancestral-parsimony.yaml)                                                                 | `ancestral`                                                           | Fitch parsimony on the graph; each segment's reconstruction equals a run on that segment's tree alone                          |
| [`100/optimize.yaml`](100/optimize.yaml), [`1564/optimize.yaml`](1564/optimize.yaml)                                                                                                             | `optimize`                                                            | Branch lengths fitted jointly to the mutations of all segments that share an edge                                              |
| [`100/timetree-confidence.yaml`](100/timetree-confidence.yaml), [`500/timetree-confidence.yaml`](500/timetree-confidence.yaml), [`1564/timetree-confidence.yaml`](1564/timetree-confidence.yaml) | `timetree`                                                            | Dated graph of the four segments with 90% date intervals                                                                       |
| [`100/timetree-ha-na.yaml`](100/timetree-ha-na.yaml), [`1564/timetree-ha-na.yaml`](1564/timetree-ha-na.yaml)                                                                                     | `timetree`                                                            | Dated graph of HA and NA, the segments of the public Nextstrain builds                                                         |
| [`100/pipeline-nextstrain.yaml`](100/pipeline-nextstrain.yaml), [`500/pipeline-nextstrain.yaml`](500/pipeline-nextstrain.yaml), [`1564/pipeline-nextstrain.yaml`](1564/pipeline-nextstrain.yaml) | `optimize` -> `timetree` -> `ancestral` -> `mugration` -> `mugration` | Nextstrain seasonal-flu H3N2 build on the graph: joint dating, ancestral sequences, and the region on the HA and NA time trees |

A run on a graph writes the files that hold sequence positions or one tree once per segment, numbered in the order of `segment_trees`: `ancestral.1.sequences.fasta` is HA, `timetree.2.nexus` is the dated NA tree, `timetree.3.auspice.json` the Auspice view of PB1. The node dates table `timetree.dates.tsv` covers every node of the graph. The pipeline reads the per-segment trees of a step by path (`{{ steps.optimize.output_all }}/optimize.1.nexus`), because a reference to an output that a step writes once per segment (`{{ steps.optimize.outputs.nexus }}`) resolves to several files and is an error.

On a graph of segment trees, `timetree` runs at a fixed clock rate, keeps the roots of the graph and its polytomies, and runs without a coalescent prior, a clock filter, rerooting and a rate uncertainty.

## Clock rate

The examples date the four segments at 0.0024 substitutions per site per year, and HA with NA at 0.0033.

Each segment of a graph gets a relative rate, estimated from the data and scaled to a mean of 1 over all alignment columns, so a branch length of the graph is in substitutions per site averaged over all sites, and the clock rate is the rate averaged over all sites. The examples take the per-segment H3N2 rates of the Nextstrain seasonal-flu workflow ([nextstrain/seasonal-flu@9ec17de, `workflow/snakemake_rules/core.smk`, lines 190-203](https://github.com/nextstrain/seasonal-flu/blob/9ec17ded9d9ae662097afccb61c087a453277014/workflow/snakemake_rules/core.smk#L190-L203), estimated from 12-year builds) and weight them by the alignment widths:

| Segment   | Columns | Nextstrain rate |      Columns x rate |
| --------- | ------: | --------------: | ------------------: |
| HA        |    1737 |         0.00382 |               6.635 |
| NA        |    1436 |         0.00267 |               3.834 |
| PB1       |    2341 |         0.00139 |               3.254 |
| PB2       |    2341 |         0.00218 |               5.103 |
| All four  |    7855 |                 | 18.827, mean 0.0024 |
| HA and NA |    3173 |                 | 10.469, mean 0.0033 |

The relative rates that the run estimates on the full set are 1.17 (HA), 1.09 (NA), 0.89 (PB1) and 0.93 (PB2). They follow the total lengths of the input trees, 1.64, 1.48, 1.25 and 1.24 substitutions per site, and differ from the ratios of the Nextstrain rates (1.59, 1.11, 0.58, 0.91). The mean rate, which sets the dates, does not depend on them. A root-to-tip regression of each input tree on its own root (`treetime clock --keep-root`) gives 0.0044, 0.0024, 0.0020 and 0.0025 for the four segments, a site-weighted mean of 0.0027.

## Clock outliers

The Nextstrain workflow removes tips whose date deviates from their inferred time by a z-score above 4 before it dates a tree, and its `augur refine` removes tips beyond 4 interquartile ranges of the root-to-tip regression. Segment graphs have no clock filter, so the datasets keep all tips. With the clock filter at 4 interquartile ranges, the root-to-tip regression of the full set on the input roots flags in HA the reference strains A/California/7/2004, A/NewYork/55/2004 and A/Wellington/1/2004 and the sample A/Malaysia/2312169/2023, and no tip in NA, PB1 and PB2. The 2004 strains are references that the Nextstrain builds always keep.

## Source

[`data/flu/treeknit/`](../../treeknit/) holds the output of TreeKnit (treeknit-rs 0.1.0, option `--better-MCCs`, gamma 2, 2 rounds, seed 1; see `parameters.json` and `log.txt.gz`) for the HA, NA, PB1 and PB2 trees of 1564 viruses that have all four segments, from a build in the layout of the Nextstrain seasonal-flu workflow named `h3n2_europe_6y`. The configuration of that build is in no repository. The sequences are open INSDC (GenBank) data, through the GenSpectrum Loculus mirror (`GS_` accessions in the `accession` column of `data/flu/treeknit/metadata.tsv.gz`).

TreeKnit resolves the polytomies of each tree against the other trees and finds, for every pair of segments, the maximum compatible clades (MCCs): sets of leaves whose subtree both trees share. Its trees name internal nodes per tree only.

## Build

[`build`](build) writes the three datasets from `data/flu/treeknit/`; the output is the same in every run. From the repository root, in the Python environment of the repository (Biopython and zstd):

```bash
pixi run python -I data/flu/h3n2/segments/build
```

1. **Node names.** [`test_scripts/treeknit_named_trees.py`](../../../../test_scripts/treeknit_named_trees.py) names the internal nodes of the four resolved trees so that an ancestor that several segments share has one name in each of their trees. For each segment pair, it matches the node of one tree with the node of the other tree that is the most recent common ancestor of the same MCC leaves, joins the matches of all six pairs into node identities, and splits the identities that hold two nodes of one tree (10 identities) or that close a cycle in the merged graph (3 identities) back into nodes of one tree. On the full set, 508 identities are shared by all four segments, 92 by three, 99 by two, and 518 internal nodes belong to one segment only. Method and limits: [`kb/reports/2026-11-15-treeknit-output-to-named-segment-trees.md`](../../../../kb/reports/2026-11-15-treeknit-output-to-named-segment-trees.md)
2. **Leaves of the subsets.** Each region gets a share of the leaves in proportion to its size (largest remainder, at least one leaf per region). Within a region, the leaves are sorted by date and taken at equal intervals, the oldest first, as [`data/subsample`](../../../subsample) does
3. **Pruning.** The named trees of the full set are pruned to the chosen leaves; nodes left without leaves are removed, and a node left with one child is merged into the edge to its child, with the summed length, as tree inference would give it. The names stay, so the nodes that the segments share stay shared: a subset of the leaves of an MCC is still compatible across the segment pair
4. **Files.** The alignments keep the sequences of the tree leaves; the metadata keeps the columns listed above

The merge of the full set inserts 179 nodes into segments that pass through them, and leaves 9 edges whose other segments disagree on the path and 29 edges whose path holds a node of the segment elsewhere (`--verbosity info` prints these counts).

## Acknowledgements

The sequences were collected, sequenced and shared by many laboratories. The `accession_ha`, `accession_na`, `accession_pb1` and `accession_pb2` columns of each `metadata.tsv` give the INSDC accession of each segment, and the `authors` and `institution` columns name the submitters. We gratefully acknowledge their work. INSDC places no restrictions on the use or redistribution of the data; some submitters may claim patent or other rights on their sequences.
