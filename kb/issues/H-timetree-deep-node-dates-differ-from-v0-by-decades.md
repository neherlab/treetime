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
- **norovirus all/genome**: the root moves by at most 0.2 years (1885.2, 1885.2, and 1885.4 at 300, 1200, and 4800 points; the run at 4800 points needs `--max-grid-points 8000000` and takes 2.5 hours). The two nodes below it move by 2.4 and 2.0 years (234 tips, 1907.0, 1904.6, 1902.6) and by 1.0 and 0.8 years (233 tips, 1936.4, 1935.4, 1934.6), so their error falls by only 1.2 times per quartering of the spacing, and no order of convergence and no Richardson estimate follow. The other clades of the larger root side among the 10 oldest nodes move by at most 1.3 years over the three grids. On the side of 34 tips, three clades agree within 0.4 years between 300 and 1200 points but move younger by 13.4 (34 tips, 1952.5 to 1965.6), 7.6 (32 tips), and 5.0 years (28 tips) at 4800 points, with 90% regions that narrow from about 2 years to 0.1 to 0.6 years

> [!IMPORTANT]
> **Investigation required.** The jump of the clades of the 34-tip root side at 4800 points is not a convergence pattern: the two coarser grids agree within 0.4 years, and the regions narrow from 1.3 to 2.2 years to 0.1 to 0.6 years. It is not shown whether the jump comes from the finer grid or from the raised `--max-grid-points`, which the coarser runs leave at its default. Evidence to collect: the run at 4800 points with the default limit raised only as far as it needs, the run at 1200 points with `--max-grid-points 8000000`, and the posteriors of the 34-tip clade in the three runs

### Coalescent prior

The coalescent prior with the actual child count is an approved v1 difference that acts on deep nodes ([kb/decisions/coalescent-total-lh-actual-multiplicity.md](../decisions/coalescent-total-lh-actual-multiplicity.md)). Without the prior in both engines, the v0-v1 root difference changes as follows:

- **lassa l**: 89.7 to 52.6 years, so the prior causes 41%. It moves v1 by 34.2 years and v0 by 2.9 years
- **rubella E1**: 16.3 to 9.7 years, so the prior causes 40%. It moves v0 by 9.5 years and v1 by 2.9 years. Without the prior, the seven clades of 1,661 to 1,672 tips among the 10 oldest nodes agree within 2.0 years, against 3.0 to 5.5 years with the prior (v0 replay); the clade of 1,674 tips differs by 4.2 years, against 11.6 years
- **norovirus all/genome**: 50.9 to 54.4 years, so the prior causes none of the difference

### Topology and branch lengths

- **lassa l**: the root splits the shared tips into 535 and 217 in both engines, and the substitution lengths (`mutation_length`) of the two root branches and the mean root-to-tip divergence of each side agree to five digits (0.04346 and 0.03245; 0.47055 and 0.30674 substitutions per site). The 52.6 years without the prior come from the time inference on the same tree and the same branch lengths
- **rubella E1**: v1 keeps the root of its input tree, whose root branch leads to one tip that the v0 clock filter removes; below it, the MRCA of the 3,023 shared tips splits them into 2,990 and 33. The v0 midpoint root splits them into 1,674 and 1,349. The deepest dates compare different splits
- **norovirus all/genome**: the root splits the tips into 234 and 34 in v1 and into 236 and 32 in v0, with root branches of 0.032 and 0.264 substitutions per site in v1 and 0.059 and 0.008 in v0. The two engines choose different roots (see [M-timetree-reroot-schedule-diverges-from-v0.md](M-timetree-reroot-schedule-diverges-from-v0.md)). Shared clades among the 10 oldest nodes differ by up to 20 years (19 tips, 1951.1 in v1 against 1971.0 in v0)

### Simulation with a known truth on lassa l

One replicate of lassa l with a known truth under a strict clock does not reproduce the difference: the two engines date the root 0.24 years apart, against 52.6 years on the real data without the prior.

- **Truth**: the topology and node dates of the v1 run without the prior (root 1539.83, root split 217 and 535), because v1 keeps the shared rooted tree as it is, while v0 resolves its polytomies
- **Data**: on each branch of time $\Delta t$, the number of substitution events is Poisson with mean $L \mu \Delta t$ ($L = 7279$ sites, $\mu = 0.00078$); each event moves a uniform site to one of the three other bases (Jukes-Cantor). The root sequence follows the base composition of the real alignment. The input tree is the true topology with branch lengths $k / L$, the event count per site that the input-mode branch likelihoods of both engines model. The metadata keeps the real date strings, so the same tips have ambiguous or missing dates
- **Runs**: the arguments of the round on the true rooted tree with the true clock rate and without the prior: v1 `--keep-root` without `--coalescent`, v0 `augur refine --keep-root --coalescent 0`. The longest input branch is 0.259 substitutions per site, so both engines use input-mode branch likelihoods, as on the real tree (longest branch 0.241)
- **Result, seed 1**: v1 1544.10 [1537.79, 1550.40], v0 1544.34 [1537.87, 1550.50]; both 90% regions contain the true date

The real input tree breaks the strict clock at the root: the mean root-to-tip divergence is 0.514 substitutions per site on the side of 535 tips and 0.339 on the side of 217 tips, a difference of about 224 years at the fixed rate, while the median tip dates of the two sides differ by 3.4 years (2013.7 and 2017.1). The simulated tree has 0.361 and 0.370. So the 52.6 years arise only on data that conflict with the strict clock, and a strict-clock simulation cannot show which engine dates lassa l correctly. The replicates of the simulation contract (100 unselected replicates and the coverage of the 90% root region) were not run, because the regime check failed.

### Chains of zero-length branches

On a caterpillar with a closed-form root age, both engines move the root younger, v0 by up to 1.29 years at 200 tips and v1 by 0.78 years at its default grid ([H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md](H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md)). This error has the same sign as the production differences (v0 younger than v1), but it is two orders of magnitude smaller, and the pre-loop collapse removes most such chains from the production trees.

## Attribution

- **lassa l**: the coalescent prior causes 37 of 90 years. The other 53 years arise in the time inference with the same topology and the same branch lengths. They need the conflict of the two root sides with the strict clock: on strict-clock data with the same topology and tip dates, the engines agree within 0.24 years. Their mechanism is not located
- **rubella E1**: the coalescent prior causes 7 of 16 years. The remaining 10 years lie at the MRCA, whose split differs between the engines, while the large clades below it agree within 2 years without the prior
- **norovirus all/genome**: the coalescent prior causes none of the 51 years, and the engines choose different roots

> [!IMPORTANT]
> **Investigation required.** The measurements do not show which engine dates the deep nodes correctly, and the 53 years of lassa l without the prior have no located cause. The evidence to collect:
>
> - **lassa l, clock conflict at the root**: how each engine weighs the two root sides, whose divergence differs by about 224 years at the fixed rate. A simulation with a known truth must reproduce this conflict (for example a rate that differs between the two sides, or a relaxed clock) before its replicates can show which engine is correct; on strict-clock data the engines agree
> - **lassa l, branch-length likelihood**: whether both engines use input-mode branch likelihoods on this tree (the root-to-tip divergence is about 0.47 substitutions per site, so many branches are longer than 0.05, where v0 uses its Gaussian cost), and the root date of v1 with the branch likelihoods of v0 evaluated on the v0 grid ([M-timetree-input-branch-likelihood-grid-diverges-from-v0.md](M-timetree-input-branch-likelihood-grid-diverges-from-v0.md))
> - **norovirus all/genome and rubella E1, rooting**: the root dates of both engines on the same rooted tree (v1 with `--keep-root` on the v0 output topology)
> - **soft-end lines**: whether the deep posteriors lie on the straight lines at soft ends of convolutions, whose slope carries FFT round-off and moves deep dates by tenths of a year ([M-distribution-convolution-soft-end-line-slope-carries-fft-round-off.md](M-distribution-convolution-soft-end-line-slope-carries-fft-round-off.md)); this error is too small to cause decades, but it bounds the precision of every comparison above

## Reproduction

The production cases of the round, with their input trees, alignments, and metadata, are in the pipelandia data packages; the v1 arguments are those of the cases above. For the grid runs, change `GRID_POINTS` and remove the compile-time check that it stays below `MaxGridPoints::MIN`. For v0 without the prior, rerun the `augur refine` command line of the build with `--coalescent 0`.

## Related

- [kb/decisions/distribution-convolution-tails-by-exponential-tilting.md](../decisions/distribution-convolution-tails-by-exponential-tilting.md): the exact tails at hard ends, which changed the v1 dates of these builds
- [kb/v0-errata/timetree-fft-convolution-rectangle-rule.md](../v0-errata/timetree-fft-convolution-rectangle-rule.md): the first-order error of v0 convolutions at hard ends
