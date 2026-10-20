# Time distributions reach far into the past at a moderate fixed clock rate

`treetime timetree` on `data/mpox/clade-ii/1000` with a fixed clock rate of `1e-3` substitutions per site per year, about ten times a realistic mpox rate, stops with an error of the grid point limit ([kb/decisions/distribution-grid-point-limit.md](../decisions/distribution-grid-point-limit.md)):

```
Time inference needs a grid of 21360606634291804 points, more than the limit of 1000000 (--max-grid-points). The fixed --clock-rate 0.001 may be far from the rate the data supports, which makes branch-time distributions very narrow. Check --clock-rate, or raise --max-grid-points when enough memory is available.
When sending the time message backward along edge 324
A grid over [-1.65682e12, 2022.62] with spacing 7.75644e-5 needs 21360606634291804 points, more than the limit of 1000000
```

The error message blames the clock rate, but the cause is the support of a node distribution, which reaches about $1.7 \cdot 10^{12}$ years into the past. A rate ten times too high does not explain such a support. Without the limit, the same run was terminated by the operating system (SIGKILL) after about 260 s at about 55 GB of resident memory. The run finishes at `1e-4` (peak about 0.9 GB, largest grid 115 156 points) and at `5e-4` (peak about 0.75 GB, largest grid 103 072 points).

## Reproduction

```bash
./dev/docker/run just r treetime timetree --tree=data/mpox/clade-ii/1000/tree.nwk \
  --metadata=data/mpox/clade-ii/1000/metadata.tsv --alignment=data/mpox/clade-ii/1000/aln.fasta.xz \
  --clock-rate=1e-3 --max-iter=1 --output-all=tmp/mass-window-explosion
```

## Evidence

Temporary logging of every grid above 10 000 points in a debug build, before the grid point limit existed, showed that `fn resample_to_mass_window()` ([packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L97-L121](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L97-L121)) computed point counts between 2.2e15 and 2.1e16 for 15 distributions, before the first convolution of the backward pass. The function takes the smaller of `width / (grid_points - 1)` and the spacing of its input (`normalized.dx()`), so these counts mean an input spacing many orders of magnitude finer than the mass window. The first calls come from `rewindow_to_mass()` at the end of `compute_branch_length_distribution()` ([packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L63](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L63)), which builds each branch-time distribution from a fixed-size branch-length grid divided by `clock_rate * gamma`.

After these requests the run continued into the backward pass, with convolution grids of up to about 37 000 points, before it was terminated.

> [!IMPORTANT]
> **Investigation required.** Collect:
>
> - The input spacing and mass window of one exploding branch distribution, and the branch length, `one_mutation`, and `gamma` of its edge, to find why the spacing collapses at `1e-3` but not at `5e-4`
> - Whether those requests allocated memory or failed, given that `rewindow_to_mass()` propagates errors with `?` and the run still continued; and which allocation reached 55 GB
> - The v0 result for the same input with `--time-marginal always`, on the decompressed alignment

## Behavior with the grid point limit

The grid point limit turns the termination into an error within about 10 s, at about 435 MB. With the limit in place, the first exceeded grid is the resampled operand of a backward-message convolution ("When sending the time message backward along edge"), whose input node distribution reaches far into the past. The `1e-3` row is measured on the current code, the other rows on earlier code:

| Clock rate | Range of the operand (years) | Spacing (years) | Points needed |
| --- | --- | --- | --- |
| `7e-4` | [-6.98e11, 2024.29] | 1.13e-4 | 6.2e15 |
| `1e-3` | [-1.66e12, 2022.62] | 7.76e-5 | 2.1e16 |
| `2e-3` | [-3.85e6, 2023.04] | 5.10e-5 | 7.6e10 |
| `0.5` | [-1.22e9, 2022.37] | 2.08e-7 | 5.9e15 |

The point counts match the windows of 2.2e15 to 2.1e16 points that the logging found, but the wide support belongs to node distributions that enter the convolution, which suggests that the support of a distribution grows across the backward pass, not only the resolution of a branch distribution.

> [!IMPORTANT]
> **Investigation required.** Find the operation that first gives a node distribution a support reaching more than a few years before the root date, at `7e-4` on `data/mpox/clade-ii/1000` with `--max-iter=1`.

## Related

- [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md): in the forward pass, hard sides are never trimmed by mass and every posterior keeps the root's lower edge
