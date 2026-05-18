# Switch skyline optimizer from Nelder-Mead to SLSQP

v1 uses Nelder-Mead (simplex) via `argmin` for skyline Tc optimization. v0 uses SLSQP via `scipy.optimize.minimize`. Both minimize the same objective: negative coalescent log-likelihood plus a stiffness penalty on `log(Tc)` differences.

- v1: `optimize_skyline()` (`#optimize_skyline`) using `NelderMead` at
  [`packages/treetime/src/coalescent/skyline.rs#L68-L126`](../../packages/treetime/src/coalescent/skyline.rs#L68-L126)
- v0: `optimize_skyline()` (`#optimize_skyline`) using `method='SLSQP'` at
  [`packages/legacy/treetime/treetime/merger_models.py#L281-L318`](../../packages/legacy/treetime/treetime/merger_models.py#L281-L318)

SLSQP is a gradient-based constrained optimizer that can exploit cost function smoothness. Nelder-Mead is gradient-free and tolerates noisy or non-smooth landscapes but converges slower and can settle in different local optima.

The skyline cost function is smooth (piecewise-linear lineage counts, quadratic stiffness penalty), making SLSQP a natural fit. Nelder-Mead was chosen for v1 because `argmin` does not include an SLSQP implementation and the `argmin` ecosystem was preferred over binding to a C/Fortran optimizer.

For datasets where the cost surface has multiple local optima, the two optimizers can converge to different solutions. For well-separated optima (typical of real phylogenetic datasets with clear population size changes), both should find the same global minimum.

## Related issues

- Source: [M-timetree-skyline-nelder-mead-optimizer.md](../issues/M-timetree-skyline-nelder-mead-optimizer.md) -- delete after full resolution
- [H-timetree-skyline-objective-mismatch.md](../issues/H-timetree-skyline-objective-mismatch.md): the skyline cost function itself diverges from the per-edge Kingman cost, independent of the optimizer choice
