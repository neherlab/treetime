# Timetree branch-length factor is mass-windowed at construction

`fn compute_branch_length_distribution()` evaluates the branch-length likelihood on a pilot grid, attaches its tail laws, and then resamples it onto a mass-bounded window with `fn rewindow_to_mass()` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L79](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L79). The stored branch-length distribution is therefore the re-windowed grid, not the pilot grid ([M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md) describes the pilot extent).

The branch-length factor has a closed form, $L(t) \propto (\mu t)^{n} e^{-\mu t}$ times the Poisson indel term, with a hard bound at $t = 0$, its peak at the maximum-likelihood branch length, and a soft log-linear tail past the peak. A grid derived from that form needs no mass step. The design of the timetree distributions assigns mass sizing to two consumers only: the output grid of a convolution across an edge, and the mass statistics of the posteriors (the 90% region and the survival functions). The construction of the branch-length factor is the one other place that mass-sizes a grid.

> [!IMPORTANT]
> **Decision required.** Keep the mass window at construction, or build the branch-length grid from the closed form of the factor without a mass step.
>
> - Keep the window: the current behavior. One mass rule bounds every stored message, and the window cuts the pilot grid where it is wider than the mass of the likelihood, but a pilot grid that is coarse against the peak can blur the peak before the window measures mass ([M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md))
> - Closed-form grid: the extent and spacing follow from $n$, $\mu$, and the indel term, so no pilot evaluation is resampled, and the peak lands on the grid by construction. It changes the branch-length grids of every timetree run, so the smoke comparison shows changed node dates, and the scientific change needs approval and an entry in `kb/decisions/`
>
> Evidence to collect for either option: the node dates and the log marginal likelihood of the timetree smoke datasets with both constructions, and the grid point counts per branch.

## Validation

- Timetree smoke comparison of both constructions, with the date changes in units of the node's 90% width
- `test_gm_runner_marginal_dense` against v0, whose gap [M-timetree-branch-grid-uniform-resolution.md](M-timetree-branch-grid-uniform-resolution.md) records
