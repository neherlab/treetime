# Per-edge optimization: Newton step clamping in reparameterized spaces

## Problem

Newton's method computes a step $\delta = \ell'/\ell''$. The step is clamped to prevent overshooting. In $t$-space, `clamp(step, -1.0, bl)` ensures $t$ does not decrease by more than 1.0 subs/site and does not go negative.

The current `newton_sqrt_inner` uses `clamp(step, -1.0, s)` -- it carries the $t$-space lower bound $-1.0$ without translating it to $s$-space. At $s = 0.3$ ($t = 0.09$), the untranslated bound allows $t$ to increase by roughly 0.4, not the intended 1.0. At other values of $s$ the discrepancy is larger or smaller. The code works by coincidence because Newton steps near the optimum are small, but the bound is incorrect and could cause divergence on extreme inputs (very short branches with many indels where the first few steps are large).

The new `newton_log_inner` also needs correct step bounds in $\ln(t)$ space.

## Correct bounds by space

### $t$-space (existing, correct -- no change)

- Upper bound $t$: prevents $t_{\text{new}} = t - \delta_t < 0$
- Lower bound $-1.0$: limits $t$-space increase to 1.0 subs/site per iteration

### $\sqrt{t}$-space

The optimization variable is $s = \sqrt{t}$. The Newton step $\delta_s$ must prevent $s < 0$ and limit the $t$-space change.

- **Upper bound** $s$: prevents $s_{\text{new}} = s - \delta_s < 0$

- **Lower bound** (the $s$-step for $t$ increasing by 1.0 subs/site):

  We want $(s - \delta_s)^2 - s^2 \leq 1$. Expanding: $\delta_s^2 - 2s\delta_s - 1 \leq 0$. The critical (most negative) $\delta_s$ is:

  $$\delta_s = s - \sqrt{s^2 + 1}$$

  So the lower bound is $-(\sqrt{s^2 + 1} - s)$, always negative (allowing $s$ to increase). This must be recomputed each iteration since $s$ changes.

### $\ln(t)$-space

The optimization variable is $u = \ln(t)$.

- **Upper bound** $u - u_{\min}$: prevents $u_{\text{new}} = u - \delta_u < u_{\min}$

- **Lower bound** (the $u$-step for $t$ increasing by 1.0 subs/site):

  We want $e^{u - \delta_u} - e^u \leq 1$, i.e., $t(e^{-\delta_u} - 1) \leq 1$. So $e^{-\delta_u} \leq 1 + 1/t$, giving:

  $$\delta_u \geq -\ln(1 + 1/t)$$

  The lower bound is $-\ln(1 + 1/t)$, always negative. Since $t = e^u$ changes each iteration, recompute per iteration.

## Approach

Fix `newton_sqrt_inner`: replace the untranslated `-1.0` lower bound with the derived $-(\sqrt{s^2+1} - s)$ expression. The upper bound $s$ is already correct.

Apply the $\ln(t)$-space bounds in `newton_log_inner`: lower bound $-\ln(1 + 1/t)$, upper bound $u - u_{\min}$.

No change to `newton_inner` ($t$-space bounds are correct).

## Verification

`./dev/docker/run ./dev/dev t` -- all tests pass. The fix tightens bounds (more conservative clamping), so existing tests that pass with the loose bounds will pass with correct bounds. The fix prevents potential divergence on inputs where the first few Newton steps are large.

## Dependencies

- Depends on: [M-optimize-method-scaffolding](M-optimize-method-scaffolding.md) -- file split puts Newton code in its own module
- Depended on by: [M-optimize-method-dispatch](M-optimize-method-dispatch.md)
- Related: [M-optimize-method-newton-log](M-optimize-method-newton-log.md) -- `newton_log_inner` uses the log-space bounds defined here

## Cross-references

- Current `newton_sqrt_inner`: `optimize_unified.rs:404` (moves to Newton module after [M-optimize-method-scaffolding](M-optimize-method-scaffolding.md))
- Current `newton_inner`: `optimize_unified.rs:359`
