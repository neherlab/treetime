# Per-edge optimization method selection: scaffolding

## Problem

v1's Newton-Raphson optimizer converges to the wrong optimum on indel-bearing edges. The Poisson indel Hessian ($-k/t^2$) dominates the substitution Hessian on short branches, making Newton steps tiny and causing the step-size convergence criterion to fire before the gradient reaches zero. The optimizer reports the indel-only MLE ($k/\mu$) instead of the true combined (substitution + indel) optimum.

v0 avoids this by using Brent's method (derivative-free, immune to Hessian conditioning) in $\sqrt{t}$ space. The fix is to provide multiple optimization methods, defaulting to one that matches v0.

## Design: why 6 methods

The problem has two orthogonal axes: algorithm (Newton vs Brent) and parameterization ($t$, $\sqrt{t}$, $\ln(t)$). A third axis (convergence criterion: step-size vs step-size + derivative-magnitude) was evaluated and excluded -- derivative-magnitude makes Newton crawl with many tiny steps instead of fixing the root cause (bad conditioning). Reparameterization is the correct fix.

All 12 combinations (2 algorithms x 3 parameterizations x 2 convergence criteria) were evaluated:

- **Newton-t**: baseline, matches RAxML-NG/IQ-TREE. Known Hessian dominance limitation on indel-bearing edges.
- **Newton-sqrt**: reduces indel singularity from $O(1/t^2)$ to $O(1/t)$. Matches v0's parameterization. Residual dominance on extreme cases ($t < 0.001$, $k > 10$).
- **Newton-log**: eliminates indel singularity entirely ($\ell''_{\text{indel}} = -\mu t$, bounded). Natural relative tolerance. Best conditioning.
- **Brent-t**: derivative-free, so conditioning-invariant. But steep $\ln(\mu t)$ curvature near $t=0$ degrades parabolic interpolation to golden section. Included for completeness (2 algorithms x 3 parameterizations = 6 methods), even though Brent-sqrt dominates it for convergence speed.
- **Brent-sqrt**: matches v0 exactly (same algorithm, same parameterization). Smooth objective. Good parabolic interpolation. Default.
- **Brent-log**: smoothest objective. Best parabolic interpolation. Requires finite lower bound in log-space.

All 6 step-size + derivative-magnitude combinations were excluded because: (a) in $t$-space it treats the symptom not the cause, (b) in $\sqrt{t}$-space marginal benefit, (c) in $\ln(t)$-space guards against a problem that does not exist. Step-size-only convergence is used for all Newton variants.

| CLI value     | Algorithm | Space      | Default? |
| :------------ | :-------- | :--------- | :------- |
| `brent`       | Brent     | $t$        |          |
| `brent-sqrt`  | Brent     | $\sqrt{t}$ | Yes      |
| `brent-log`   | Brent     | $\ln(t)$   |          |
| `newton`      | Newton    | $t$        |          |
| `newton-sqrt` | Newton    | $\sqrt{t}$ |          |
| `newton-log`  | Newton    | $\ln(t)$   |          |

Default is `brent-sqrt` because it matches v0 exactly, enabling golden master comparison against v0 reference outputs.

## This issue: two preparation steps

### Step 1: Split `optimize_unified.rs`

`optimize_unified.rs` is 743 lines with 24 functions. Adding more functions (Brent parameterized, Newton log, step clamping) would push it past 800. The clock command uses a similar split: `method_brent.rs`, `method_golden_section.rs`, `method_grid.rs` at `packages/treetime/src/commands/clock/find_best_root/`.

Proposed split into 3 files:

- Newton functions (chain rules, inner loops, tolerance, constants) into a new Newton module
- Brent functions (inner loop, cost function struct, constants) into a new Brent module
- Shared infrastructure stays (metrics, contributions, evaluation, grid search, zero-branch logic, `run_optimize_mixed`, `initial_guess_mixed`)

Pure refactor, no behavior change. All functions that move become `pub(crate)`.

### Step 2: Rewrite `BranchOptMethod` enum

Current state: 3 variants (`NewtonSqrt` default, `Newton`, `Brent`-in-t-space) in a separate `optimize_method.rs` file.

Target: 6 variants with `BrentSqrt` as default. The single-enum file is unnecessary -- the clock command's `OptimizationMethod` lives in `params.rs` alongside related structs. Move the enum to `args.rs` where `InitialGuessMode` and `TreetimeOptimizeArgs` already live.

The `--opt-method` CLI help text needs updating to describe all 6 methods.

New variants (`BrentSqrt`, `BrentLog`, `NewtonLog`) can use `todo!()` in the dispatch match until the Brent, Newton-log, and step-clamping issues implement them. No test calls these variants yet.

## Verification

`./dev/docker/run ./dev/dev t` after each step -- all tests pass, no behavior change.

## Dependencies

None. This is the first step. All other method selection issues depend on this.

## Cross-references

- Next: [M-optimize-method-brent-parameterized](M-optimize-method-brent-parameterized.md), [M-optimize-method-newton-log](M-optimize-method-newton-log.md), [M-optimize-method-step-clamping](M-optimize-method-step-clamping.md)
- Clock command split pattern: `packages/treetime/src/commands/clock/find_best_root/`
- Current code: `packages/treetime/src/commands/optimize/optimize_unified.rs`, `optimize_method.rs`, `args.rs`
