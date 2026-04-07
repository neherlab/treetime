# Newton premature convergence when indel Hessian dominates substitution curvature

The Newton optimizer in `run_optimize_mixed()` converges on step size, not derivative magnitude. When the Poisson indel second derivative is much larger in magnitude than the substitution second derivative, Newton steps become tiny and the loop exits near the indel-only optimum ($k / \mu$) instead of the true combined (substitution + indel) optimum.

## Mechanism

The Newton step is $\delta = \ell'(t) / \ell''(t)$. The convergence criterion at [packages/treetime/src/commands/optimize/optimize_unified.rs#L420](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L420) is step-size based:

$$|t_{\text{new}} - t_{\text{old}}| < \max(0.001 \cdot t_{\text{old}},\; 10^{-8})$$

The Poisson indel second derivative is $-k / t^2$ (from [packages/treetime/src/commands/optimize/optimize_indel.rs#L36](../../packages/treetime/src/commands/optimize/optimize_indel.rs#L36)), which grows as $O(1/t^2)$ near zero. For short branches with indels, this term dominates the total Hessian.

When the total Hessian is dominated by the indel term, the Newton step $\delta = \ell'_{\text{total}} / \ell''_{\text{total}}$ becomes small even when $\ell'_{\text{total}}$ is large. The step falls below the convergence tolerance before the derivative reaches zero.

## Observed behavior

With 4 indels on an edge (2 indels per partition, dense + sparse), the optimizer converges to $t = 0.090$:

| Component    | $\ell'(t)$ | $\ell''(t)$ |
| :----------- | ---------: | ----------: |
| Substitution |      -18.8 |         -55 |
| Indel        |        0.0 |        -493 |
| **Combined** |  **-18.8** |    **-548** |

Newton step: $-18.8 / -548 = 0.034$. Tolerance at $t = 0.09$: $\max(0.001 \times 0.09, 10^{-8}) \approx 9 \times 10^{-5}$. The step (0.034) exceeds tolerance so the loop continues, but subsequent steps shrink rapidly as the indel Hessian increases near the indel MLE. The loop converges to $t = k / \mu_{\text{indel}} = 4 / 44.4 = 0.090$ -- the Poisson MLE where the indel derivative is zero but the substitution derivative is -18.8.

The combined log-likelihood at $t = 0.072$ is -10.16, higher than -10.40 at the reported "optimum" $t = 0.090$. The optimizer did not find a local maximum of the combined objective.

## Location

- Newton loop: [packages/treetime/src/commands/optimize/optimize_unified.rs#L413-L435](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L413-L435)
- Newton tolerance: [packages/treetime/src/commands/optimize/optimize_unified.rs#L32-L39](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L32-L39) (`NEWTON_REL_TOL = 0.001`, `NEWTON_ABS_TOL = 1e-8`)
- Poisson indel Hessian: [packages/treetime/src/commands/optimize/optimize_indel.rs#L36](../../packages/treetime/src/commands/optimize/optimize_indel.rs#L36) (`-k / t^2`)

## v0 comparison

v0 uses `scipy.optimize.minimize_scalar` with Brent's method at [packages/legacy/treetime/treetime/gtr.py#L881-L891](../../packages/legacy/treetime/treetime/gtr.py#L881-L891). Brent is derivative-free and bracket-based, finding the true minimum within a bracket. v0 also reparameterizes in $\sqrt{t}$ space (the objective evaluates at $t^2$, the result is squared back), which flattens the objective near zero and improves conditioning.

The v1 Newton implementation has no v0 precedent. v0 does not use Newton's method for branch length optimization.

## Impact

Medium. On edges where indel curvature dominates (many indels relative to branch length, or short branches), the optimized branch length is biased toward the indel-only MLE. The substitution signal is effectively ignored at convergence. The bias depends on the ratio $|\ell''_{\text{indel}}| / |\ell''_{\text{sub}}|$, which is $O(k / (t^2 \cdot |\ell''_{\text{sub}}|))$.

In practice, phylogenetic datasets typically have few indels per edge (k = 0-3) and moderate branch lengths (t = 0.01-0.1), so the indel Hessian dominance is limited to short branches with multiple indels. For standard flu/ebola/zika datasets, this may not produce visible bias. For datasets with high indel rates (e.g. some bacterial genomes), the bias could be significant.

## Possible approaches

### A1. Replace Newton with Brent's method

Match v0's approach. The `argmin` crate provides `BrentOpt` (or use `BrentRoot` for the derivative-zero problem). The bracket is available: lower bound from `grid_search_branch_lengths()`, upper bound from the same. Brent is derivative-free and guaranteed to find the optimum within the bracket.

- Pro: eliminates the Hessian dominance problem entirely; matches v0's proven approach; robust to ill-conditioned Hessians
- Con: derivative-free methods converge slower than Newton when the Hessian is well-conditioned; requires bracket construction (already available from grid search path); loses the diagnostic information from derivatives

The $\sqrt{t}$ reparameterization from v0 could also be adopted to improve conditioning near zero.

### A2. Add derivative-magnitude convergence criterion

Require both step size AND derivative magnitude to be small before accepting convergence. For example: exit only when $|t_{\text{new}} - t_{\text{old}}| < \text{tol}$ AND $|\ell'_{\text{total}}(t)| < \text{deriv\_tol}$.

- Pro: minimal code change; preserves Newton's quadratic convergence when the Hessian is well-conditioned
- Con: requires choosing `deriv_tol`, which depends on the scale of the log-likelihood (varies with alignment length); may cause non-convergence if the derivative tolerance is too tight

### A3. Line search with Wolfe conditions

Replace the simple Newton step with a damped Newton step that satisfies the Wolfe conditions (sufficient decrease + curvature condition). This ensures each step actually decreases the objective.

- Pro: guaranteed decrease per step; standard optimization practice; handles ill-conditioned Hessians gracefully
- Con: more complex implementation; requires objective evaluation at trial points; the `argmin` crate provides `MoreThuenteLineSearch` but integration requires restructuring the loop

### A4. Rescale indel contribution

Scale the indel log-likelihood and derivatives by a factor that balances the Hessian magnitudes. For example, divide the indel contribution by $\max(1, |\ell''_{\text{indel}}| / |\ell''_{\text{sub}}|)$.

- Pro: keeps Newton; simple change
- Con: changes the objective function, so the optimum shifts; requires theoretical justification that the rescaled objective produces scientifically correct branch lengths; the scaling factor depends on quantities that change during iteration

## Related

- [Dense and sparse evaluation loops duplicated](L-optimize-eval-dense-sparse-duplication.md) -- resolved eval loop duplication in the same code region
- [Indel rate estimation bypassed in Never mode on all-zero-BL trees](N-optimize-indel-rate-never-mode-zero-bl.md) -- related indel rate estimation edge case
