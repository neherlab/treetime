# `simulated/reassortment/80`

A simulated ancestral reassortment graph of a genome with four segments, with 80 samples and a known truth: the true date of every node, the true clock rate and the true substitution models. Use it to try `ancestral`, `optimize` and `timetree` on graphs of segment trees (`--segment-trees`, `--segment-alignments`) and to check their results against the truth.

## Simulation

One replicate of the reassortment simulation of [validation/segments.py](../../../../validation/segments.py) (`simulate` step, model in [validation/reassortment.py](../../../../validation/reassortment.py)):

- **Graph**: coalescent with reassortment (CoalRe model) with a constant coalescent time scale of 10 years and a reassortment rate of 0.2 per lineage and year (case `high`). 2 samples at each of 40 sample times, 0.1 years apart, the most recent dated 2025.0. The graph has 27 reassortments
- **Segments**: each with its own length, substitution rate and random GTR model; sequences evolve along the segment's tree of the graph

  | Segment | Sites | Rate (substitutions per site and year) |
  | ------- | ----: | -------------------------------------: |
  | `s1`    |  1500 |                                  0.004 |
  | `s2`    |  1200 |                                  0.003 |
  | `s3`    |  1000 |                                  0.002 |
  | `s4`    |   800 |                                  0.001 |

- **True clock rate**: 0.00275556 substitutions per site and year, the mean of the segment rates weighted by segment length. The examples date the graph at this rate
- **Input trees**: the collapsed segment trees, as tree inference gives them: no nodes that a segment only passes through, and no hidden nodes (nodes with one child in every segment, such as the reassortment nodes themselves, which no tree inference recovers). Every node keeps its true name, so the trees merge into the true graph, as after an ideal matching of the nodes across the trees. Branch lengths are divergences, the true rate times the true time span times a lognormal(0, 0.3) factor per branch

## Files

| File                      | Content                                                                                                                                                       |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `tree_s{1,2,3,4}.nwk`     | Collapsed segment trees with named internal nodes; branch lengths in substitutions per site                                                                   |
| `aln_s{1,2,3,4}.fasta.xz` | Sample sequences of each segment                                                                                                                              |
| `metadata.tsv`            | Sample dates                                                                                                                                                  |
| `truth_nodes.tsv`         | Every node of the merged trees: its true numeric date, and its number of parents in the true graph without hidden nodes (2 or 3 for 20 nodes, 0 for the root) |
| `truth_model.json`        | True clock rate, reassortment rate, number of reassortments, and per segment the length, rate, equilibrium frequencies (`pi`) and exchange rates (`W`)        |
| `*.yaml`                  | Example configs; the header of each file gives the command line that runs it                                                                                  |

## Expected accuracy

[kb/reports/2026-11-15-segment-validation-timetree.md](../../../../kb/reports/2026-11-15-segment-validation-timetree.md) measures `timetree --segment-trees` on 3 replicates of this case: on collapsed trees, the dates of internal nodes have a bias below 0.14 years and an RMSE of 0.25 to 0.64 years, and their 90% regions cover 85 to 89% of the true dates. On this replicate, `timetree.yaml` gives an RMSE of about 0.3 years over the 101 internal nodes.

## Regenerate

From `validation/`, in the pixi environment of the repository:

```sh
pixi run python segments.py simulate --reps 1 --seed 2026 --data ../tmp/segments
pixi run python segments.py export --data ../tmp/segments --case high --rep 0 --out ../data/simulated/reassortment/80
```

The simulation code is that of commit `d9d3c6ba5`; `export` copies the collapsed trees and the alignments of the replicate and writes `metadata.tsv` and the truth files. The seed and the replicate fix every random draw, so the commands rewrite the same files.
