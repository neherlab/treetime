# Deep-node dates differ from v0 by decades on three production builds

On three Nextstrain production builds, v0 and v1 date the oldest nodes 16 to 90 years apart, with v1 older in every case. The measurements below exclude two of the suspected causes. The discretization error of v1 is not a cause, because the v1 dates are converged in the time grid. The error along chains of zero-length branches is not a cause, because it moves dates younger and by at most about one year. The coalescent prior causes about 40% of the difference on two builds and none on the third. Which engine dates these nodes correctly is open.

## Cases

The v1 runs use the `timetree` arguments of the pipelandia round (`--max-iter 2 --clock-filter 0 --confidence --branch-length-mode auto --divergence-units mutations-per-site --seed 1`, the fixed clock rate of each workflow, one thread). v0 is `augur refine` 34.3.0 with TreeTime 0.12.1 in the Nextstrain image of the round, with the command line of the workflow:

- **lassa l**: 752 tips, clock rate 0.00078, `--coalescent opt`; v1 `--keep-root` on the midpoint-rooted tree that augur writes, v0 `--root mid_point`
- **rubella E1**: 3,046 tips in v1 and 3,023 in v0 (the v0 clock filter `--clock-filter-iqd 4` removes 23), clock rate 0.0006, `--coalescent opt`; v1 `--keep-root --clock-filter 4.0 --prune-filtered-tips`, v0 `--root mid_point --stochastic-resolve`
- **norovirus all/genome** (example build of the CI profile): 268 tips, clock rate 0.00328; v1 `--coalescent 0.01 --reroot least-squares`, v0 `--root best --stochastic-resolve` without `--coalescent`, which augur runs with $T_c = 0.01$

Nodes are matched by their tip sets restricted to the tips that all runs share: each restricted set maps to its most recent common ancestor in each tree. Dates are the numeric dates that the runs report, with the 90% max posterior region. v1 reports the posterior median at the root and the posterior peak at other nodes.

## Evidence

### Root dates

v1 at `85b5eb96a` (exact convolution tails at hard ends). "v0 replay" reruns the v0 command line of the round and reproduces the frozen v0 root dates within 0.5 years and the 10 oldest nodes below them within 3.5 years. "No coalescent" removes the coalescent prior from both engines: v1 without `--coalescent`, v0 with `--coalescent 0`, which augur passes to TreeTime as $T_c = 0$, so no prior.

| Build | v0 frozen | v0 replay | v1 | v0, no coalescent | v1, no coalescent |
| --- | ---: | ---: | ---: | ---: | ---: |
| lassa l | 1595.8 [1585.7, 1608.4] | 1595.3 | 1505.6 [1491.5, 1518.7] | 1592.4 | 1539.8 |
| rubella E1, MRCA of the shared tips | 1919.9 [1917.5, 1925.0] | 1919.7 | 1903.4 [1896.2, 1910.0] | 1910.2 | 1900.5 |
| norovirus all/genome | 1936.2 [1934.5, 1938.2] | 1936.1 | 1885.2 [1881.6, 1888.8] | 1934.4 | 1880.0 |

### Grid convergence of v1

v1 with `GRID_POINTS` (`packages/treetime/src/timetree/inference/runner.rs`) at 300 (default) and 1200 points. This constant sets the points of the branch distributions, the backward and forward convolution outputs, and the coalescent factor:

- **lassa l**: the root and the 10 oldest internal nodes move by at most 0.4 years (root 1505.6 and 1505.8)
- **rubella E1**: at most 0.2 years (MRCA 1903.4 and 1903.6)
- **norovirus all/genome**: the root does not move (1885.2), but the two nodes below it move by 2.4 years (234 tips, 1907.0 to 1904.6) and 1.0 year (233 tips, 1936.4 to 1935.4), and the other 8 of the 10 oldest nodes by at most 0.8 years

### Coalescent prior

The coalescent prior with the actual child count is an approved v1 difference that acts on deep nodes ([kb/decisions/coalescent-total-lh-actual-multiplicity.md](../decisions/coalescent-total-lh-actual-multiplicity.md)). Without the prior in both engines, the v0-v1 root difference changes as follows:

- **lassa l**: 89.7 to 52.6 years, so the prior causes 41%. It moves v1 by 34.2 years and v0 by 2.9 years
- **rubella E1**: 16.3 to 9.7 years, so the prior causes 40%. It moves v0 by 9.5 years and v1 by 2.9 years. Without the prior, the seven clades of 1,661 to 1,672 tips among the 10 oldest nodes agree within 2.0 years, against 3.0 to 5.5 years with the prior (v0 replay); the clade of 1,674 tips differs by 4.2 years, against 11.6 years
- **norovirus all/genome**: 50.9 to 54.4 years, so the prior causes none of the difference

### Topology and branch lengths

- **lassa l**: the root splits the shared tips into 535 and 217 in both engines, and the substitution lengths (`mutation_length`) of the two root branches and the mean root-to-tip divergence of each side agree to five digits (0.04346 and 0.03245; 0.47055 and 0.30674 substitutions per site). The 52.6 years without the prior come from the time inference on the same tree and the same branch lengths
- **rubella E1**: v1 keeps the root of its input tree, whose root branch leads to one tip that the v0 clock filter removes; below it, the MRCA of the 3,023 shared tips splits them into 2,990 and 33. The v0 midpoint root splits them into 1,674 and 1,349. The deepest dates compare different splits
- **norovirus all/genome**: the root splits the tips into 234 and 34 in v1 and into 236 and 32 in v0, with root branches of 0.032 and 0.264 substitutions per site in v1 and 0.059 and 0.008 in v0. The two engines choose different roots (see [M-timetree-reroot-schedule-diverges-from-v0.md](M-timetree-reroot-schedule-diverges-from-v0.md)). Shared clades among the 10 oldest nodes differ by up to 20 years (19 tips, 1951.1 in v1 against 1971.0 in v0)

### Chains of zero-length branches

On a caterpillar with a closed-form root age, both engines move the root younger, v0 by up to 1.29 years at 200 tips and v1 by 0.78 years at its default grid ([H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md](H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md)). This error has the same sign as the production differences (v0 younger than v1), but it is two orders of magnitude smaller, and the pre-loop collapse removes most such chains from the production trees.

## Attribution

- **lassa l**: the coalescent prior causes 37 of 90 years. The other 53 years arise in the time inference with the same topology and the same branch lengths. Their cause is not located
- **rubella E1**: the coalescent prior causes 7 of 16 years. The remaining 10 years lie at the MRCA, whose split differs between the engines, while the large clades below it agree within 2 years without the prior
- **norovirus all/genome**: the coalescent prior causes none of the 51 years, and the engines choose different roots

> [!IMPORTANT]
> **Investigation required.** The measurements do not show which engine dates the deep nodes correctly, and the 53 years of lassa l without the prior have no located cause. The evidence to collect:
>
> - **Simulation with a known truth**: 100 unselected replicates with recorded seeds, with the true topology and clock rate given to both engines, and the exact Clopper-Pearson 95% interval of the coverage of the 90% root date interval of each engine
> - **lassa l, branch-length likelihood**: whether both engines use input-mode branch likelihoods on this tree (the root-to-tip divergence is about 0.47 substitutions per site, so many branches are longer than 0.05, where v0 uses its Gaussian cost), and the root date of v1 with the branch likelihoods of v0 evaluated on the v0 grid ([M-timetree-input-branch-likelihood-grid-diverges-from-v0.md](M-timetree-input-branch-likelihood-grid-diverges-from-v0.md))
> - **norovirus all/genome and rubella E1, rooting**: the root dates of both engines on the same rooted tree (v1 with `--keep-root` on the v0 output topology)
> - **soft-end lines**: whether the deep posteriors lie on the straight lines at soft ends of convolutions, whose slope carries FFT round-off and moves deep dates by tenths of a year ([M-distribution-convolution-soft-end-line-slope-carries-fft-round-off.md](M-distribution-convolution-soft-end-line-slope-carries-fft-round-off.md)); this error is too small to cause decades, but it bounds the precision of every comparison above

## Reproduction

The production cases of the round, with their input trees, alignments, and metadata, are in the pipelandia data packages; the v1 arguments are those of the cases above. For the grid runs, change `GRID_POINTS` and remove the compile-time check that it stays below `MaxGridPoints::MIN`. For v0 without the prior, rerun the `augur refine` command line of the build with `--coalescent 0`.

## Related

- [kb/decisions/distribution-convolution-tails-by-exponential-tilting.md](../decisions/distribution-convolution-tails-by-exponential-tilting.md): the exact tails at hard ends, which changed the v1 dates of these builds
- [kb/v0-errata/timetree-fft-convolution-rectangle-rule.md](../v0-errata/timetree-fft-convolution-rectangle-rule.md): the first-order error of v0 convolutions at hard ends
