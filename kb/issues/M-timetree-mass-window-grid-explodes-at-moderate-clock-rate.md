# Mass-window resampling requests astronomically large grids at a moderate fixed clock rate

`treetime timetree` on `data/mpox/clade-ii/1000` with a fixed clock rate of `1e-3` substitutions per site per year, about ten times a realistic mpox rate, is terminated by the operating system (SIGKILL) after about 260 s at about 55 GB of resident memory. The same run finishes at `1e-4` (peak about 0.9 GB, largest grid 115 156 points) and at `5e-4` (peak about 0.75 GB, largest grid 103 072 points).

## Reproduction

```bash
./dev/docker/run just r treetime timetree --tree=data/mpox/clade-ii/1000/tree.nwk \
  --metadata=data/mpox/clade-ii/1000/metadata.tsv --alignment=data/mpox/clade-ii/1000/aln.fasta.xz \
  --clock-rate=1e-3 --max-iter=1 --output-all=tmp/mass-window-explosion
```

## Evidence

Temporary logging of every grid above 10 000 points in a debug build showed that `fn resample_to_mass_window()` ([packages/treetime-distribution/src/distribution_ops/mass_domain.rs](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs)) computed point counts between 2.2e15 and 2.1e16 for 15 distributions, before the first convolution of the backward pass. The function takes the smaller of `width / (grid_points - 1)` and the spacing of its input (`normalized.dx()`), so these counts mean an input spacing many orders of magnitude finer than the mass window. The first calls come from `rewindow_to_mass()` at the end of `compute_branch_length_distribution()` ([packages/treetime/src/timetree/inference/branch_length_likelihood.rs](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs)), which builds each branch-time distribution from a fixed-size branch-length grid divided by `clock_rate * gamma`.

After these requests the run continued into the backward pass, with convolution grids of up to about 37 000 points, before it was terminated.

> [!IMPORTANT]
> **Investigation required.** Collect:
>
> - The input spacing and mass window of one exploding branch distribution, and the branch length, `one_mutation`, and `gamma` of its edge, to find why the spacing collapses at `1e-3` but not at `5e-4`
> - Whether those requests allocated memory or failed, given that `rewindow_to_mass()` propagates errors with `?` and the run still continued; and which allocation reached 55 GB
> - The v0 result for the same input with `--time-marginal always`, on the decompressed alignment

## Related

- [M-timetree-marginal-dense-mpox-slow.md](M-timetree-marginal-dense-mpox-slow.md): grid growth from spacing ratios in the convolution on the same organism
