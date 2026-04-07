# Per-edge branch length optimizer converges to wrong optimum on indel-bearing edges

The Newton-Raphson optimizer in `run_optimize_mixed()` converges prematurely when the Poisson indel second derivative dominates the substitution second derivative. The convergence criterion checks only step size ($|t_{\text{new}} - t_{\text{old}}| < \text{tol}$), not derivative magnitude. When the combined Hessian is dominated by the indel term ($-k/t^2$), Newton steps shrink below tolerance while the first derivative of the substitution component is still large. The optimizer reports the indel-only MLE ($k/\mu$) instead of the true combined (substitution + indel) optimum.

End-to-end test coverage for indel-aware optimization is incomplete: existing tests assert only positivity and finiteness of the result, not local optimality of the combined log-likelihood.

## Per-edge likelihood function

The per-site conditional likelihood for one edge, using eigendecomposition $Q = V \, \text{diag}(\lambda) \, V^{-1}$, is:

$$L_i(t) = \sum_c k_{ic} \exp(\lambda_c t)$$

where $k_{ic}$ are precomputed coefficients depending on endpoint messages (independent of $t$) and $\lambda_c$ are the substitution model eigenvalues. The total substitution log-likelihood sums over sites:

$$\ell_{\text{sub}}(t) = \sum_i \ln L_i(t)$$

The first and second derivatives follow by multiplying $k_{ic}$ by $\lambda_c$ and $\lambda_c^2$ respectively, reusing the same $\exp(\lambda_c t)$ values:

$$\ell'_{\text{sub}}(t) = \sum_i \frac{\sum_c k_{ic} \lambda_c \exp(\lambda_c t)}{\sum_c k_{ic} \exp(\lambda_c t)}$$

$$\ell''_{\text{sub}}(t) = \sum_i \left[ \frac{\sum_c k_{ic} \lambda_c^2 \exp(\lambda_c t)}{L_i(t)} - \left(\ell'_{i}(t)\right)^2 \right]$$

Dense implementation: [packages/treetime/src/commands/optimize/optimize_dense.rs#L67-L75](../../packages/treetime/src/commands/optimize/optimize_dense.rs#L67-L75). Sparse implementation: [packages/treetime/src/commands/optimize/optimize_sparse.rs#L46-L119](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L46-L119).

## Poisson indel model

v1 models indels on an edge as a Poisson process with rate $\mu$ per unit branch length. Given $k$ observed indel events on a branch of length $t$:

$$\ell_{\text{indel}}(t) = k \ln(\mu t) - \mu t - \ln(k!)$$

$$\ell'_{\text{indel}}(t) = \frac{k}{t} - \mu$$

$$\ell''_{\text{indel}}(t) = -\frac{k}{t^2}$$

The MLE is at $t_{\text{indel}} = k / \mu$. When $k = 0$: $\ell(t) = -\mu t$, $\ell'(t) = -\mu$, $\ell''(t) = 0$.

Implementation: `poisson_indel_log_lh()` at [packages/treetime/src/commands/optimize/optimize_indel.rs#L20-L39](../../packages/treetime/src/commands/optimize/optimize_indel.rs#L20-L39).

The indel Hessian $-k/t^2$ has a $1/t^2$ singularity at the origin. This is intrinsic to the Poisson model: any count model where $P(k > 0 | t) \to 0$ as $t \to 0$ must have $\ell(t) \to -\infty$ near zero, producing divergent curvature.

## Hessian dominance mechanism

The total log-likelihood is:

$$\ell(t) = \ell_{\text{sub}}(t) + \ell_{\text{indel}}(t)$$

The substitution Hessian scales with alignment length: $\ell''_{\text{sub}} \sim O(L)$. The indel Hessian scales as $O(k/t^2)$. For short branches with indels, $|\ell''_{\text{indel}}| \gg |\ell''_{\text{sub}}|$.

The Newton step is:

$$\delta = \frac{\ell'(t)}{\ell''(t)} = \frac{\ell'_{\text{sub}} + \ell'_{\text{indel}}}{\ell''_{\text{sub}} + \ell''_{\text{indel}}}$$

When $|\ell''_{\text{indel}}| \gg |\ell''_{\text{sub}}|$, the denominator is dominated by the indel term. Near the indel MLE where $\ell'_{\text{indel}} \approx 0$, the numerator is $\approx \ell'_{\text{sub}}$ (which can be large), but the denominator is $\approx \ell''_{\text{indel}}$ (very large in magnitude). The step $\delta$ becomes small even though $\ell'_{\text{sub}}$ is large.

The convergence criterion at [packages/treetime/src/commands/optimize/optimize_unified.rs#L420](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L420) is:

$$|t_{\text{new}} - t_{\text{old}}| < \max(0.001 \cdot t_{\text{old}},\; 10^{-8})$$

Small Newton steps satisfy this tolerance before the combined gradient reaches zero. The optimizer exits at the indel MLE, not the combined optimum.

### Numerical example

With 4 indels on an edge (2 indels per partition, dense + sparse), $\mu = 44.4$, the optimizer converges to $t = 0.090$:

| Component    | $\ell'(t)$ | $\ell''(t)$ |
| :----------- | ---------: | ----------: |
| Substitution |    $-18.8$ |       $-55$ |
| Indel        |      $0.0$ |      $-493$ |
| **Combined** |    $-18.8$ |      $-548$ |

Newton step: $-18.8 / -548 = 0.034$. Subsequent steps shrink rapidly as the indel Hessian increases. The loop converges to $t = k/\mu = 4/44.4 = 0.090$ where the indel derivative is zero but the substitution derivative is $-18.8$.

The combined log-likelihood at $t = 0.072$ is $-10.16$, higher than $-10.40$ at the reported "optimum" $t = 0.090$.

## Root cause analysis

The core issue is a conditioning problem in the combined objective, not a defect in Newton's method:

- Substitution Hessian: $O(L)$ where $L$ is alignment length
- Indel Hessian: $O(k/t^2)$, diverging at $t = 0$

Newton's step size is controlled by the total Hessian. When one component dominates, Newton preferentially tracks that component's optimum. The step-size convergence criterion then fires before the other component's gradient is zeroed.

RAxML-NG's convergence criterion checks $|\ell'(t)| < 10^{-4}$ (derivative magnitude), not just step size. IQ-TREE uses a minimum branch length criterion. v1 is the only tool among RAxML-NG, IQ-TREE, PhyML, and v0 that checks step size alone.

## Unimodality of the 1D phylogenetic likelihood

<a id="cite-1a"></a>[Dinh and Matsen 2017](https://doi.org/10.1214/16-AAP1240) [[1](#ref-1)] proved:

- **JC69, F81, binary symmetric models**: the 1D likelihood has at most one stationary point. Any optimizer converges to the global MLE from any starting point.
- **K2P, HKY, GTR**: explicit counterexamples with two local maxima exist.

For TreeTime's common use case (JC69 on viral data), Newton-Raphson is provably safe for the substitution component alone. The addition of the indel component (a concave function) preserves unimodality when the substitution component is unimodal. For GTR models, the existing grid search fallback (triggered when $\ell''(t) \geq 0$) mitigates multimodality risk.

## Fix: $\sqrt{t}$ reparameterization

Define $s = \sqrt{t}$, optimize $\ell(s^2)$ over $s$. By the chain rule:

$$\frac{d\ell}{ds} = 2s \cdot \frac{d\ell}{dt}$$

$$\frac{d^2\ell}{ds^2} = 4s^2 \cdot \frac{d^2\ell}{dt^2} + 2 \cdot \frac{d\ell}{dt}$$

For the indel component:

$$\frac{d^2\ell_{\text{indel}}}{ds^2} = 4t \cdot \left(-\frac{k}{t^2}\right) + 2\left(\frac{k}{t} - \mu\right) = -\frac{2k}{t} - 2\mu$$

The singularity drops from $O(1/t^2)$ to $O(1/t)$. At $t = 0.09$, $k = 4$, $\mu = 44.4$:

| Space      | $\ell''_{\text{indel}}$ | Ratio to $\ell''_{\text{sub}}$ |
| :--------- | ----------------------: | -----------------------------: |
| $t$        |                  $-493$ |                            9.0 |
| $\sqrt{t}$ |                  $-178$ |                  $\approx 3.2$ |

Newton in $s$-space no longer has its step magnitude controlled primarily by the indel curvature. The reparameterization does not change the optimum location -- it changes the coordinate system, improving conditioning of the combined Hessian.

v0 `optimal_t_compressed()` at [packages/legacy/treetime/treetime/gtr.py#L816-L920](../../packages/legacy/treetime/treetime/gtr.py#L816-L920) uses Brent's method in $\sqrt{t}$ space. v0 does not have an indel model for branch length optimization -- its reparameterization addresses general near-zero conditioning, not the specific indel Hessian dominance.

## Implementation: CLI-selectable methods

Three per-edge methods via `--opt-method`, following the clock command's `--branch-split-method` pattern:

- **`newton-sqrt`** (default): Newton-Raphson in $\sqrt{t}$ space. Transforms $t$-space derivatives via chain rule. Reduces indel Hessian dominance from $O(1/t^2)$ to $O(1/t)$. Convergence criterion in $s$-space: $|s_{\text{new}} - s_{\text{old}}| < \max(0.001 \cdot s,\; 10^{-8})$.

- **`newton`**: Newton-Raphson in $t$ space. Subject to Hessian dominance on short branches with indels.

- **`brent`**: Brent's method via `argmin::BrentOpt`. Derivative-free, bracket-based. Bracket: $[\text{min\_bl},\; \max(1.5t + 1/L,\; 0.5)]$. Immune to Hessian conditioning. Convergence independent of curvature ratios.

All methods share: grid search fallback for non-concave regions, zero-branch short-circuit for unimodal models.

## Newton-Raphson convergence properties

Newton-Raphson converges quadratically near a non-degenerate stationary point: each iteration doubles the number of correct digits. For the per-edge phylogenetic likelihood, typical convergence is 3-8 iterations. This is the reason RAxML-NG, IQ-TREE, and PAML use Newton-Raphson as the primary per-branch optimizer (<a id="cite-2a"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[2](#ref-2)]; <a id="cite-3a"></a>[Stamatakis 2014](https://doi.org/10.1093/bioinformatics/btu033) [[3](#ref-3)]; <a id="cite-4a"></a>[Nguyen et al. 2015](https://doi.org/10.1093/molbev/msu300) [[4](#ref-4)]).

The eigendecomposition-based coefficient caching makes derivatives nearly free: $\ell'(t)$ and $\ell''(t)$ reuse the same $\exp(\lambda_c t)$ values as $\ell(t)$.

## Brent's method

<a id="cite-5a"></a>[Brent 1973](https://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf) [[5](#ref-5)] combines golden section search with successive parabolic interpolation. Convergence order is approximately 1.325 (tribonacci constant). Brent is invariant to conditioning: it finds the optimum within the bracket regardless of curvature ratios. v0 uses Brent; RAxML-NG and IQ-TREE keep it as a fallback for Newton.

## Test gaps to address

1. **End-to-end optimality verification**: tests must verify the reported branch length is at a local maximum of the combined log-likelihood (C1: log-likelihood at optimum exceeds neighbors at multiple delta scales), not just assert positivity and finiteness.

2. **Stationarity**: the implied Newton step at the optimum ($|\ell'/\ell''|$) should be smaller than the convergence tolerance (C2).

3. **Cross-method agreement**: NewtonSqrt and Brent should achieve similar combined log-likelihood values on the same input (C3).

4. **Brent bracket validity**: the log-likelihood at the reported optimum must exceed the log-likelihood at both bracket endpoints (C4).

5. **NewtonSqrt improvement over Newton**: on the Hessian-dominated case, NewtonSqrt must achieve equal or better log-likelihood than Newton in $t$-space (C5).

6. **Chain rule correctness**: $\sqrt{t}$-space analytical derivatives must match numerical finite differences.

7. **Zero-branch mismatch precondition**: the guard `branch_length == 0.0 && !all_sites_valid_at_zero()` may be unreachable through marginal reconstruction with JC69. If confirmed, document as defense-in-depth.

8. **Finite-difference Hessian tolerance**: the proptest second-derivative check uses `max_relative = 1e-3`, which is 1000x looser than the theoretical bound for central differences.

## Location

- Newton loop: [packages/treetime/src/commands/optimize/optimize_unified.rs#L413-L435](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L413-L435)
- Newton tolerance: [packages/treetime/src/commands/optimize/optimize_unified.rs#L32-L39](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L32-L39) (`NEWTON_REL_TOL = 0.001`, `NEWTON_ABS_TOL = 1e-8`)
- Poisson indel Hessian: [packages/treetime/src/commands/optimize/optimize_indel.rs#L36](../../packages/treetime/src/commands/optimize/optimize_indel.rs#L36) (`-k / t^2`)

## v0 comparison

v0 uses `scipy.optimize.minimize_scalar` with Brent's method in $\sqrt{t}$ space at [packages/legacy/treetime/treetime/gtr.py#L881-L891](../../packages/legacy/treetime/treetime/gtr.py#L881-L891). v0 does not have a Poisson indel model for branch length optimization. The v1 Newton implementation has no v0 precedent.

## Impact

Medium. On edges where indel curvature dominates (many indels relative to branch length, or short branches), the optimized branch length is biased toward the indel-only MLE. The bias depends on the ratio $|\ell''_{\text{indel}}| / |\ell''_{\text{sub}}|$, which is $O(k / (t^2 \cdot |\ell''_{\text{sub}}|))$. In practice, phylogenetic datasets typically have few indels per edge ($k = 0$-$3$) and moderate branch lengths ($t = 0.01$-$0.1$), so the bias is limited to short branches with multiple indels. For datasets with high indel rates (some bacterial genomes), the bias could be significant.

## Related

- [docs/port-proposals/optimize-convergence-and-method-choice.md](../../docs/port-proposals/optimize-convergence-and-method-choice.md) -- broader proposal (P1-P6); P5 (Brent as alternative) is subsumed here
- [docs/port-intentional-changes/optimize-newton-raphson-per-edge.md](../../docs/port-intentional-changes/optimize-newton-raphson-per-edge.md) -- rationale for v1 choosing Newton over Brent
- [docs/reports/optimization-methods/2-branch-length-optimization.md](../../docs/reports/optimization-methods/2-branch-length-optimization.md) -- tool comparison: Newton vs Brent across RAxML-NG, IQ-TREE, PhyML, v0, v1
- [L-optimize-eval-dense-sparse-duplication.md](L-optimize-eval-dense-sparse-duplication.md) -- enum dispatch duplication in the same code region

## Supersedes

- `M-optimize-newton-indel-hessian-dominance.md` (the production convergence bug)
- `N-optimize-indel-newton-test-gaps.md` (test coverage gaps blocked by the convergence bug)

## References

1. <a id="ref-1"></a> Dinh, Vu, and Frederick A. Matsen IV. 2017. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." _Ann. Appl. Prob._ 27(3):1250-1286. https://doi.org/10.1214/16-AAP1240 [↩](#cite-1a)
2. <a id="ref-2"></a> Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _J. Mol. Evol._ 17:368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-2a)
3. <a id="ref-3"></a> Stamatakis, Alexandros. 2014. "RAxML Version 8: A Tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies." _Bioinformatics_ 30(9):1312-1313. https://doi.org/10.1093/bioinformatics/btu033 [↩](#cite-3a)
4. <a id="ref-4"></a> Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh. 2015. "IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies." _Mol. Biol. Evol._ 32(1):268-274. https://doi.org/10.1093/molbev/msu300 [↩](#cite-4a)
5. <a id="ref-5"></a> Brent, Richard P. 1973. _Algorithms for Minimization without Derivatives._ Prentice-Hall. https://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf [↩](#cite-5a)
