# `simulated/reassortment/80`

A simulated ancestral reassortment graph of a genome with four segments, with 80 samples and a known truth: the true date of every node, the true clock rate and the true substitution models. Use it to try `ancestral`, `optimize` and `timetree` on graphs of segment trees (`--segment-trees`, `--segment-alignments`) and to check their results against the truth.

## Simulation

One replicate of the reassortment simulation of [validation/segments.py](../../../../validation/segments.py) (`simulate` step, model in [validation/reassortment.py](../../../../validation/reassortment.py)):

- **Graph**: coalescent with reassortment (CoalRe model) with a constant coalescent time scale of 10 years (case `high`). A lineage that carries more than one segment draws a reassortment at a rate of 0.2 per year, and each of its segments goes to one of two parents with probability 1/2; a draw that sends all segments to the same parent changes nothing, so the rate of actual reassortments is lower. 2 samples at each of 40 sample times, 0.1 years apart, the most recent dated 2025.0. The graph has 27 reassortments
- **Segments**: each with its own length, substitution rate and random GTR model; sequences evolve along the segment's tree of the graph

  | Segment | Sites | Rate (substitutions per site and year) |
  | ------- | ----: | -------------------------------------: |
  | `s1`    |  1500 |                                  0.004 |
  | `s2`    |  1200 |                                  0.003 |
  | `s3`    |  1000 |                                  0.002 |
  | `s4`    |   800 |                                  0.001 |

- **True clock rate**: 0.00275556 substitutions per site and year, the mean of the segment rates weighted by segment length. The examples date the graph at this rate
- **Input trees**: per segment, the true tree of the segment without the nodes that the segment only passes through, as tree inference gives it from that segment's sequences: it starts at the most recent common ancestor of the segment, and it has no hidden nodes (nodes with one child in every segment, such as the reassortment nodes themselves, which no tree inference recovers). Every node keeps its true name, as after an ideal matching of the nodes across the trees. Branch lengths are divergences, the true rate times the true time span times a lognormal(0, 0.3) factor per branch
- **Inputs and true graph**: the merged input trees differ from the true graph where a removed pass-through node cannot be placed from the other trees. In the true graph, `node_0132` (root of the `s2` tree) descends from `node_0133` (root of the `s1` and `s3` trees), and `node_0128` (root of the `s4` tree) from `node_0132`, along edges on which these segments only pass through; no input tree holds these edges, so the merged inputs have three roots against one. Merged by TreeTime (`--output-tree-dot` of `optimize`), the inputs also give 26 nodes with several parents, against 20 in the true graph

## Files

| File                      | Content                                                                                                                                                       |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `tree_s{1,2,3,4}.nwk`     | Collapsed segment trees with named internal nodes; branch lengths in substitutions per site                                                                   |
| `aln_s{1,2,3,4}.fasta.xz` | Sample sequences of each segment                                                                                                                              |
| `metadata.tsv`            | Sample dates                                                                                                                                                  |
| `truth_nodes.tsv`         | Every node of the input trees (the same nodes as the true graph without hidden nodes): `numeric_date`, its true date; `parents`, its number of parents in the true graph with the hidden nodes contracted (2 or 3 for 20 nodes, 0 for the root `node_0133`) |
| `truth_model.json`        | Seed, case and replicate of the simulation; true clock rate, reassortment draw rate, number of actual reassortments, and per segment the length, rate, equilibrium frequencies (`pi`) and exchange rates (`W`) |
| `*.yaml`                  | Example configs; the header of each file gives the command line that runs it                                                                                  |

## Expected accuracy

[kb/reports/2026-11-15-segment-validation-timetree.md](../../../../kb/reports/2026-11-15-segment-validation-timetree.md) measures `timetree --segment-trees` on 3 replicates of this case: on collapsed trees, the dates of internal nodes have a bias below 0.14 years and an RMSE of 0.25 to 0.64 years, and their 90% regions cover 85 to 89% of the true dates. On this replicate, `timetree.yaml` gives an RMSE of about 0.3 years over the 101 internal nodes.

## Regenerate

From `validation/`, in the pixi environment of the repository:

```sh
data="$(mktemp -d)"
pixi run python segments.py simulate --reps 1 --seed 2026 --data "${data}"
pixi run python segments.py export --seed 2026 --data "${data}" --case high --rep 0 --out ../data/simulated/reassortment/80
```

The seed and the replicate fix every random draw; the draws are those of the simulation code at commit `d9d3c6ba5`, so the commands rewrite the same files. `simulate` records the seed and a hash of every file it writes, and `export` stops when the seed or a file it reads differs, because `simulate` keeps a replicate that a directory already holds. `export` copies the collapsed trees and the alignments of the replicate and writes `metadata.tsv` and the truth files.
