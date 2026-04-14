# Per-edge optimization: Newton-Raphson in ln(t) space

## Problem

Newton-Raphson in $t$-space converges to the wrong optimum when the Poisson indel Hessian dominates. The $\sqrt{t}$ reparameterization (existing `newton_sqrt_inner`) reduces the singularity from $O(1/t^2)$ to $O(1/t)$ but does not eliminate it. For extreme cases ($t < 0.001$, $k > 10$), residual dominance remains.

The $\ln(t)$ reparameterization eliminates the indel singularity entirely, producing bounded indel curvature ($-\mu t$) regardless of branch length.

## Scientific background

### Per-edge likelihood

The per-site conditional likelihood for one edge, using eigendecomposition $Q = V \operatorname{diag}(\lambda) V^{-1}$:

$$L_i(t) = \sum_c k_{ic} \exp(\lambda_c t)$$

where $k_{ic}$ are precomputed coefficients (independent of $t$) and $\lambda_c$ are eigenvalues. With Poisson indel model ($k$ observed indels, rate $\mu$):

$$\ell_{\text{indel}}(t) = k \ln(\mu t) - \mu t - \ln(k!), \quad \ell''_{\text{indel}}(t) = -k/t^2$$

### Hessian dominance

Example at $t = 0.09$, $k = 4$, $\mu = 44.4$, $\ell''_{\text{sub}} \approx -55$:

| Space      | $\ell''_{\text{indel}}$ | Ratio to $\ell''_{\text{sub}}$ |
| :--------- | ----------------------: | -----------------------------: |
| $t$        |                  $-493$ |                          $9.0$ |
| $\sqrt{t}$ |                  $-178$ |                          $3.2$ |
| $\ln(t)$   |                  $-4.0$ |                         $0.07$ |

In $\ln(t)$ space the indel term is 14x smaller than the substitution term. The dominance is inverted.

### Chain rule for ln(t) space

Define $u = \ln(t)$, optimize $\ell(e^u)$:

$$d\ell/du = t \cdot d\ell/dt$$

$$d^2\ell/du^2 = t^2 \cdot d^2\ell/dt^2 + t \cdot d\ell/dt$$

The indel Hessian in $u$-space becomes $-\mu t$ (bounded, no singularity). The $1/t^2$ singularity is eliminated entirely, not just reduced.

Derivation: $dt/du = e^u = t$, so $d\ell/du = (d\ell/dt) \cdot t$. For the second derivative: $d^2\ell/du^2 = d/du(t \cdot d\ell/dt) = t \cdot d\ell/dt + t^2 \cdot d^2\ell/dt^2$.

### Existing chain_rule_sqrt for comparison

$s = \sqrt{t}$: $d\ell/ds = 2s \cdot d\ell/dt$, $d^2\ell/ds^2 = 4s^2 \cdot d^2\ell/dt^2 + 2 \cdot d\ell/dt$.

Implemented at `packages/treetime/src/commands/optimize/method_newton.rs`.

## Approach

Three changes, all in the Newton module:

### 1. Chain rule log function

Analogous to the existing `chain_rule_sqrt`. Takes $t$-space derivatives, returns $u$-space derivatives using the formulas above. The $O(1)$ per-evaluation cost is a few multiplications.

### 2. Newton log-space tolerance

The existing `newton_tolerance` computes `max(rel_tol * x, abs_tol)`. This assumes the optimization variable is non-negative. In $\ln(t)$ space, $u = \ln(t) < 0$ for all $t < 1$ (nearly all branch lengths), making the relative term negative and collapsing the tolerance to the absolute floor.

The tolerance function should be refactored to work with any parameterization. One approach: use `x.abs()` in the relative term, with per-space constants that can be tuned independently. The constants can start identical across spaces.

A tolerance of $\epsilon$ in $u$-space corresponds to approximately $\epsilon$ relative tolerance in $t$-space (since $dt/t \approx du$ for small $du$). So a tolerance of $10^{-3}$ in $u$-space is approximately 0.1% relative change in $t$.

The tolerance refactoring affects all three Newton variants (`newton_inner`, `newton_sqrt_inner`, `newton_log_inner`) since the function signature changes.

### 3. Newton log inner loop

Same structure as `newton_sqrt_inner` (use it as template) but using the log chain rule and $u = \ln(t)$ as the optimization variable.

Key differences from `newton_sqrt_inner`:

- Initial value: $u = \ln(t)$ (requires $t > 0$; the indel starting-point logic in `run_optimize_mixed` ensures this)
- Min bound: $u_{\min} = \ln(\max(\text{min\_bl}, 10^{-12}))$ (finite, unlike $\sqrt{t}$ where $s_{\min} = 0$ is valid)
- Each iteration: evaluate fresh metrics at $t = e^u$, transform via chain rule
- Step clamping: see [M-optimize-method-step-clamping](M-optimize-method-step-clamping.md) for correct bounds
- Back-transform: $t = e^{u_{\text{final}}}$
- Grid search fallback when the $u$-space Hessian is non-negative (same pattern as other Newton variants)

Convergence uses step-size only. Derivative-magnitude was evaluated and excluded: it makes Newton crawl with many tiny steps instead of fixing the root cause (bad conditioning). The $\ln(t)$ reparameterization eliminates the indel singularity, making derivative-magnitude checks unnecessary.

### Chain rule tests

Add tests mirroring the existing `chain_rule_sqrt` test pattern:

- Analytical values at known $(t, d\ell/dt, d^2\ell/dt^2)$ inputs
- Behavior at small $t$ (both derivatives approach zero because the $t$ factor suppresses them)
- Numerical first-derivative agreement (central difference of $\ell(e^u)$) parameterized by $(t, k, \mu)$
- Numerical second-derivative agreement, parameterized

Example analytical check at $t = 0.09$, $d\ell/dt = 10$, $d^2\ell/dt^2 = -100$:

- $d\ell/du = 0.09 \times 10 = 0.9$
- $d^2\ell/du^2 = 0.09^2 \times (-100) + 0.09 \times 10 = -0.81 + 0.9 = 0.09$

## Verification

`./dev/docker/run ./dev/dev t` -- chain rule tests pass immediately, `newton_log_inner` is not called by dispatch until [M-optimize-method-dispatch](M-optimize-method-dispatch.md) wires it.

## Dependencies

- Depends on: M-optimize-method-scaffolding (done) -- file split puts Newton code in its own module
- Depended on by: [M-optimize-method-dispatch](M-optimize-method-dispatch.md)
- Related: [M-optimize-method-step-clamping](M-optimize-method-step-clamping.md) -- defines the step clamping bounds used by the log-space inner loop

## Cross-references

- Existing `chain_rule_sqrt`: `packages/treetime/src/commands/optimize/method_newton.rs`
- Existing `newton_sqrt_inner`: `packages/treetime/src/commands/optimize/method_newton.rs`
- Existing `newton_tolerance`: `packages/treetime/src/commands/optimize/optimize_unified.rs`
