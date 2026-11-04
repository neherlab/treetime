# Deep-node dates differ from v0 by decades in both directions

On the time trees of the Nextstrain workflows, v0 and v1 date the oldest nodes up to 175 years apart, in both directions: with a fixed clock rate and the same root split, v1 is older by 6 to 90 years on 14 trees and younger by 5 to 108 years on 7 others (see "Root dates across the Nextstrain workflows"). The measurements below on three builds (lassa l, rubella E1, norovirus all/genome) exclude two of the suspected causes. The discretization error of v1 is not a cause of the root differences, because the v1 root dates of all three builds move by at most 0.4 years between grids of 300 and 1200 points. Below the root of norovirus, the dates of some deep clades are not converged: they move by up to 13.4 years at 4800 points (see "Grid convergence of v1"). The error along chains of zero-length branches is not a cause, because it moves dates younger and by at most about one year. The coalescent prior causes about 40% of the difference on two builds and none on the third. Which engine dates these nodes correctly is open.

## Cases

The v1 runs use these `timetree` arguments (`--max-iter 2 --clock-filter 0 --confidence --branch-length-mode auto --divergence-units mutations-per-site --seed 1`, the fixed clock rate of each workflow, one thread). v0 is `augur refine` 34.3.0 with TreeTime 0.12.1 in the Nextstrain image of the round, with the command line of the workflow:

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

The real input tree breaks the strict clock at the root: the mean root-to-tip divergence is 0.514 substitutions per site on the side of 535 tips and 0.339 on the side of 217 tips, a difference of about 224 years at the fixed rate, while the median tip dates of the two sides differ by 3.4 years (2013.7 and 2017.1). The simulated tree has 0.361 and 0.370. So in the one replicate the 52.6 years do not arise on strict-clock data, and a strict-clock simulation cannot show which engine dates lassa l correctly. The replicates of the simulation contract (100 unselected replicates and the coverage of the 90% root region) were not run, because the regime check failed.

### Root dates across the Nextstrain workflows

The comparison of TreeTime v0 and v1 in the Nextstrain workflows (pipelandia, round 3) ran v1 at `0e522d97e` and v0 (`augur refine` 34.3.0 with TreeTime 0.12.1) with the same commands on the same inputs in the same image, on the published builds and on the example data of each workflow. The table lists the trees whose root dates differ by more than 5 years; "v1 minus v0" is negative where v1 is older. Where one engine removes tips that the other keeps, the dates are those of the common ancestor of the shared tips.

| Class | Trees, v1 minus v0 in years |
| --- | --- |
| Fixed rate, same root split, v1 older | lassa l -90.2, gpc -69.6, s -60.9; yellow-fever genome -70.0; hmpv all/genome -42.0; norovirus GII.4/genome -29.4; WNV all-lineages -15.5, lineage-2 -6.3; example data: yellow-fever genome -40.9, nipah all -17.4, WNV all-lineages -13.6, dengue all/E -11.7, lassa l -9.2, lassa s -7.9 |
| Fixed rate, same root split, v1 younger | norovirus all/p22 +108.1, all/VP1 +48.0; seasonal-cov hku1 +39.5, nl63 +7.7; oropouche L +5.2; example data: seasonal-cov hku1 +38.8, norovirus all/VP1 +18.3 |
| Fixed rate, different root split | mpox clade I +174.7; norovirus all/genome +158.7, all/RdRp +60.6, GII.17/genome -8.2; dengue all/genome -18.7; rubella E1 -17.0, genome -5.3; example data: norovirus all/genome -51.0 |
| Estimated rate | measles genome/global +86.4 (v1 rate 1.76 times v0), N450/global -38.6 (0.75); mumps global +14.7 (1.20); rsv b/F/6y -10.0 (0.68), b/genome/all-time -6.1 (0.79); example data: measles N450/global -28.8 (0.83), mumps global +8.8 (1.17), avian-flu h5n1/ha/all-time +5.5 (1.10) |

- **v1 younger**: the four norovirus builds of all genotypes set `--root best`. With an older v1 root on the same split, they would follow the pattern of lassa l; instead v1 is 48 to 159 years younger. On all/p22 the v1 root moved from 1699.8 at the v1 commit of round 1 to 1660.9
- **seasonal-cov hku1**: the engines keep different tips, because v1 removes 18 tips that the v0 filter of `augur refine` keeps ([M-clock-filter-kept-tips-differ-from-augur-refine-which-filters-without-the-clock-rate.md](M-clock-filter-kept-tips-differ-from-augur-refine-which-filters-without-the-clock-rate.md))
- **Estimated rate**: the root difference follows the rate ratio: v1 is older where its rate is lower, and younger where it is higher

### Chains of zero-length branches

On a caterpillar with a closed-form root age, both engines move the root younger, v0 by up to 1.29 years at 200 tips and v1 by 0.78 years at its default grid ([H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md](H-timetree-backward-messages-lose-accuracy-along-chains-of-zero-length-branches.md)). This error has the same sign as the production differences (v0 younger than v1), but it is two orders of magnitude smaller, and the pre-loop collapse removes most such chains from the production trees.

## Attribution

- **lassa l**: the coalescent prior causes 37 of 90 years. The other 53 years arise in the time inference with the same topology and the same branch lengths. In one replicate (seed 1) of a strict-clock simulation with the same topology and tip dates, the engines agree within 0.24 years, so the difference may need the conflict of the two root sides with the strict clock; the 100 replicates of the simulation contract were not run. The mechanism of the 53 years is not located
- **rubella E1**: the coalescent prior causes 7 of 16 years. The remaining 10 years lie at the MRCA, whose split differs between the engines, while the large clades below it agree within 2 years without the prior
- **norovirus all/genome**: the coalescent prior causes none of the 51 years, and the engines choose different roots

> [!IMPORTANT]
> **Investigation required.** The measurements do not show which engine dates the deep nodes correctly, and the 53 years of lassa l without the prior have no located cause. v0 is no ground truth, so each comparison answers one of three questions with its own evidence: do both engines solve the same problem, does each engine compute the answer of its model, and is the answer close to the true dates. The evidence to collect, the first item before the others:
>
> - **Equal inputs**: for each compared build, the tips, root, branch lengths, coalescent prior, branch-length mode, and polytomy resolution that each engine uses, traced from the source of augur and of the v1 call. A date difference between runs with different inputs or models does not show an error of the time inference; known input differences are the clock filter of `augur refine` and the rooting of `--root best`
> - **Independent numerical oracle**: the posterior of the root date on small trees by direct numerical integration of the model, in code that neither engine shares, because the model defines the correct answer
> - **Simulation with a known truth**: a regime check first, which shows that the simulated data have the conflicts of the real data (rates that differ between the two root sides, a relaxed clock), then 100 seeded replicates per case, with the bias, the spread, and the coverage of the 90% root region of both engines
> - **Tip-date cross-validation**: on the real data, hide the dates of dated tips, infer them, and compare the inferred dates with the hidden dates, which measures the accuracy of each engine without a model truth
> - **v1 younger on the norovirus builds of all genotypes and mpox clade I**: whether far tails of the v0 convolutions, which v0 crops and bridges with straight lines, cause these differences; this explanation is an expectation, not a measurement
> - **lassa l, clock conflict at the root**: how each engine weighs the two root sides, whose divergence differs by about 224 years at the fixed rate. A simulation with a known truth must reproduce this conflict (for example a rate that differs between the two sides, or a relaxed clock) before its replicates can show which engine is correct; in one strict-clock replicate (seed 1) the engines agree
> - **lassa l, branch-length likelihood**: whether both engines use input-mode branch likelihoods on this tree (the root-to-tip divergence is about 0.47 substitutions per site, so many branches are longer than 0.05, where v0 uses its Gaussian cost), and the root date of v1 with the branch likelihoods of v0 evaluated on the v0 grid ([M-timetree-input-branch-likelihood-grid-diverges-from-v0.md](M-timetree-input-branch-likelihood-grid-diverges-from-v0.md))
> - **norovirus all/genome and rubella E1, rooting**: the root dates of both engines on the same rooted tree (v1 with `--keep-root` on the v0 output topology)
> - **Stability of v1**: root dates of v1 that move by years between v1 commits ([M-timetree-root-dates-move-years-between-v1-commits-without-a-traced-cause.md](M-timetree-root-dates-move-years-between-v1-commits-without-a-traced-cause.md)) limit every comparison of these builds
> - **soft-end lines**: whether the deep posteriors lie on the straight lines at soft ends of convolutions, whose slope carries FFT round-off and moves deep dates by tenths of a year ([M-distribution-convolution-soft-end-line-slope-carries-fft-round-off.md](M-distribution-convolution-soft-end-line-slope-carries-fft-round-off.md)); this error is too small to cause decades, but it bounds the precision of every comparison above

## Reproduction

The inputs are the outputs of Nextstrain workflow runs, with their input trees, alignments, and metadata:

- **lassa l**: the `phylogenetic-build-output` artifact 12836482720 of a scheduled GitHub Actions build of `nextstrain/lassa` at commit `8bd66e384a7c458332e6d3ca07c4619efa349aa0`
- **rubella E1**: the `phylogenetic-build-output` artifact 11115633017 of `nextstrain/rubella` at commit `5d94165bbbb69a7fefd89d29bc7eb1df9ac2e5d5`
- **norovirus all/genome**: the example data and `phylogenetic/build-configs/ci/config.yaml` of `nextstrain/norovirus` as of 2026-10-28, run with the Nextstrain image of the round

GitHub deletes Actions artifacts after their retention period, so the two artifacts can expire. The v1 arguments are those of the cases above. For the grid runs, change `GRID_POINTS` and remove the compile-time check that it stays below `MaxGridPoints::MIN`. For v0 without the prior, rerun the `augur refine` command line of the build with `--coalescent 0`.

## Related

- [kb/decisions/distribution-convolution-tails-by-exponential-tilting.md](../decisions/distribution-convolution-tails-by-exponential-tilting.md): the exact tails at hard ends, which changed the v1 dates of these builds
- [kb/v0-errata/timetree-fft-convolution-rectangle-rule.md](../v0-errata/timetree-fft-convolution-rectangle-rule.md): the first-order error of v0 convolutions at hard ends
- [M-clock-filter-kept-tips-differ-from-augur-refine-which-filters-without-the-clock-rate.md](M-clock-filter-kept-tips-differ-from-augur-refine-which-filters-without-the-clock-rate.md): the tips that `augur refine` keeps differ from the tips of the v1 filter
- [M-timetree-root-dates-move-years-between-v1-commits-without-a-traced-cause.md](M-timetree-root-dates-move-years-between-v1-commits-without-a-traced-cause.md): v1 root dates that move between v1 commits
