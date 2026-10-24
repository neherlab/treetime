# Timetree backward message holds NaN on lassa/L/500

The first backward pass of the time inference stops on `lassa/L/500` for two families of settings:

```
   0: When sending the time message backward along edge 712
   1: Cannot re-peak a distribution on [1391.9880057889266, 1988.8991965630894]: its values contain NaN
```

A distribution that holds NaN is a numerical defect, whatever the input: the run must either finish or stop with a message that names the input problem.

## Reproduction

The input root kept, with a fixed clock rate and a coalescent prior:

```bash
treetime timetree --tree=data/lassa/L/500/tree.nwk --metadata=data/lassa/L/500/metadata.tsv \
  --metadata-id-columns=accession --alignment=data/lassa/L/500/aln.fasta.xz \
  --keep-root --clock-rate=0.00078 --coalescent=opt --seed=42 --output-all=<dir>
```

The `timetree/lassa/L/500/keep-root-fixed-clock` row of `dev/smoke.toml` runs this command, and `data/lassa/L/500/pipeline-nextstrain.yaml` (midpoint root, fixed rate, clock filter off) fails the same way.

| Settings | Result |
| --- | --- |
| `--keep-root --clock-rate=0.00078` | finishes |
| `--keep-root --clock-rate=0.00078 --confidence` | finishes |
| `--keep-root --clock-rate=0.00078 --no-indels` | finishes |
| `--keep-root --clock-rate=0.00078 --coalescent=opt` | NaN on edge 712 |
| `--keep-root --clock-rate=0.00078 --clock-filter=0` | NaN on edge 752 |
| `--reroot=midpoint` | finishes |
| `--reroot=midpoint --clock-rate=0.00078` | finishes |
| `--reroot=midpoint --clock-rate=0.00078 --coalescent=opt` | finishes |
| `--reroot=midpoint --clock-rate=0.00078 --confidence` | finishes |
| `--reroot=midpoint --clock-rate=0.00078 --clock-filter=0` | NaN on edge 826 |
| `--reroot=midpoint --coalescent=opt --confidence --clock-filter=0` (estimated rate) | NaN on edge 826 |

The failures are the same for every seed tried (1 to 5 and 42), and with `--keep-polytomies`, so they happen before polytomy resolution. With the default rerooting and the default clock filter, the runs finish.

## Observations

- `--branch-length-mode auto` chooses `input`, because the longest input branch is 0.218 substitutions per site
- The input root is the root that IQ-TREE wrote, not a root chosen for the temporal signal. With the kept root and the fixed rate, the clock filter marks 56 outliers, many with apparent dates near the year 1100. With the clock filter off, every such tip enters the fit and the node times
- With the kept root, the coalescent prior is seeded with one segment over [1184.9, 2022.7] and $T_c = 2351$
- The failure comes in the first backward pass, while the branch distributions are built from the input lengths

> [!IMPORTANT]
> **Investigation required.** Find the operation that first produces NaN on edges 712, 752, and 826 (input-mode branch likelihood, the product with the coalescent factor, or the convolution), and the input of that operation. Then decide whether the cause is a numerical defect on a valid input or an input that the run must reject with a message.

## Related issues

- [M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md](M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md): `lassa/L/20` stops at the grid point limit in input mode
- [N-distribution-mixed-nan-policy-undecided.md](N-distribution-mixed-nan-policy-undecided.md): no common rule for NaN in distribution operations
