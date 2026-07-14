# Skyline objective and optimizer diverge from v0

Let $T_c$ denote the coalescent population-size time scale. v1 uses Nelder-Mead via `argmin` for skyline $T_c$ optimization, while v0 uses SLSQP via `scipy.optimize.minimize`. Their objective functions, bounds, and initialization also differ, so changing only the optimizer does not establish parity.

- v1: `optimize_skyline()` (`#optimize_skyline`) using `NelderMead` at
  [`packages/treetime/src/coalescent/skyline.rs#L68-L126`](../../packages/treetime/src/coalescent/skyline.rs#L68-L126)
- v0: `optimize_skyline()` (`#optimize_skyline`) using `method='SLSQP'` at
  [`packages/legacy/treetime/treetime/merger_models.py#L281-L318`](../../packages/legacy/treetime/treetime/merger_models.py#L281-L318)

V1 clamps `log_tc` before evaluating both likelihood and stiffness. V0 evaluates stiffness on the raw vector and applies its lower-bound penalty separately. The lower penalty, initialization, and handling of the non-smooth bound also differ. Nelder-Mead and SLSQP can then converge to different local optima on already-different surfaces.

## Decision required

- Decide whether exact v0 objective behavior is required for clamping, stiffness, lower-bound penalty, and initialization.
- Once the objective contract is fixed, select an optimizer that reproduces v0 outputs within the numerical contract.
- Define golden masters covering active bounds, multiple starting points, and non-smooth or locally multimodal surfaces.

No implementation ticket is ready until the objective and bound contracts are approved.
