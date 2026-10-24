# Timetree backward message holds NaN with a kept root and a fixed rate on lassa/L/500

The first backward pass of the time inference stops on `lassa/L/500` when the input root is kept and the clock rate is fixed, together with a coalescent prior or with the clock filter switched off:

```
   0: When sending the time message backward along edge 712
   1: Cannot re-peak a distribution on [1391.9880057889266, 1988.8991965630894]: its values contain NaN
```

A distribution that holds NaN is a numerical defect, whatever the input: the run must either finish or stop with a message that names the input problem.

## Reproduction

```bash
treetime timetree --tree=data/lassa/L/500/tree.nwk --metadata=data/lassa/L/500/metadata.tsv \
  --metadata-id-columns=accession --alignment=data/lassa/L/500/aln.fasta.xz \
  --keep-root --clock-rate=0.00078 --coalescent=opt --seed=42 --output-all=<dir>
```

The `timetree/lassa/L/500/keep-root-fixed-clock` row of `dev/smoke.toml` runs this command.

| Settings added to `--keep-root --clock-rate=0.00078` | Result |
| --- | --- |
| none | finishes |
| `--confidence` | finishes |
| `--no-indels` | finishes |
| `--coalescent=opt` | NaN on edge 712 |
| `--clock-filter=0` | NaN on edge 752 |
| `--coalescent=opt --confidence --clock-filter=0 --no-indels` (the Nextstrain lassa settings) | NaN on edge 712 |

The failure is the same for every seed tried (1 to 5 and 42), and with `--keep-polytomies`, so it happens before polytomy resolution. Without `--keep-root`, or with an estimated rate, the same runs finish.

## Observations

- `--branch-length-mode auto` chooses `input`, because the longest input branch is 0.218 substitutions per site
- The input root is the root that IQ-TREE wrote, not a root chosen for the temporal signal. With the fixed rate, the clock filter marks 56 outliers, many with apparent dates near the year 1100
- The coalescent prior is seeded with one segment over [1184.9, 2022.7] and $T_c = 2351$
- The failure comes in the first backward pass, while the branch distributions are built from the input lengths

> [!IMPORTANT]
> **Investigation required.** Find the operation that first produces NaN on edge 712 (input-mode branch likelihood, the product with the coalescent factor, or the convolution), and the input of that operation. Then decide whether the cause is a numerical defect on a valid input or an input that the run must reject with a message.

## Related issues

- [M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md](M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md): `lassa/L/20` stops at the grid point limit in input mode
- [N-distribution-mixed-nan-policy-undecided.md](N-distribution-mixed-nan-policy-undecided.md): no common rule for NaN in distribution operations
