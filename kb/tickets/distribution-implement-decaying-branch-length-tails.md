# Implement decaying far tails for branch-length distributions

Branch-length time distributions are discretized on a hard-truncated grid and then evaluated beyond that grid as a flat `Constant` tail. The flat tail assigns a non-decaying probability floor to arbitrarily long branches, which biases deep internal-node dates by up to tens of years. Replace the flat far tail with a log-linear (exponential-in-probability) decay that continues the branch-length likelihood's own slope, so the far region reflects the physical decay of long branches instead of a flat floor.

This supersedes the flat-`Constant` rationale recorded in [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md) for the past direction and resolves [kb/proposals/exponential-branch-length-tails.md](../proposals/exponential-branch-length-tails.md).

## Problem

### Where the tail comes from

The branch-length time distribution is built on a grid that runs from a small positive floor to `max_bl = min(5 * center, 5.0)` substitutions/site, where `center` is the current maximum-likelihood branch length [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L62-L81](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L62-L81). The grid is a hard truncation of the long-branch side; nothing represents branches longer than `max_bl`.

The grid is converted to a time axis and returned as a `Distribution::Function` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L52-L58](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L52-L58).

In the backward pass, the branch-length distribution is negated and convolved with the child's time distribution to produce the parent message [packages/treetime/src/timetree/inference/backward_pass.rs#L124-L127](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime/src/timetree/inference/backward_pass.rs#L124-L127). The convolution result grid is exactly the sum-support of its operand grids [packages/treetime-distribution/src/distribution_ops/convolve.rs#L179-L190](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime-distribution/src/distribution_ops/convolve.rs#L179-L190) -- no tail is added during convolution. The far-past reach of the parent message is therefore bounded by the truncated `max_bl`, and beyond that boundary the message is evaluated as a flat `Constant` value [packages/treetime-grid/src/grid_fn.rs#L388-L399](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime-grid/src/grid_fn.rs#L388-L399).

The flat `Constant` left tail on the backward message is consumed by multiplication when sibling messages are combined [packages/treetime/src/timetree/inference/backward_pass.rs#L80-L84](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime/src/timetree/inference/backward_pass.rs#L80-L84): the `Constant` tail lets a message be evaluated at its boundary value in the region where a sibling's grid extends further into the past [packages/treetime-distribution/src/distribution_ops/multiply.rs#L319-L349](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime-distribution/src/distribution_ops/multiply.rs#L319-L349).

### Why the flat tail is wrong

Under a Poisson substitution model, the likelihood of a branch as a function of its length $t$ with $k$ observed substitutions is

$$\mathcal{L}(t) \;\propto\; e^{-\mu t}\,(\mu t)^{k},$$

where $t$ is branch length (substitutions/site or, after rescaling, time), $\mu$ is the per-site substitution rate, and $k$ is the number of observed substitutions on the branch. This is the branch contribution of the pruning likelihood <a id="cite-1a"></a>[Felsenstein 1981](https://doi.org/10.1007/bf01734359) [[1](#ref-1)]. For $t$ beyond the mode its logarithm is asymptotically linear, $\ln \mathcal{L}(t) \sim -\mu t$, so the density decays exponentially. A flat tail replaces this decay with a constant floor $\mathcal{L}(t_{\max})$ that never decreases, over-weighting arbitrarily long branches and, after negation and convolution, over-weighting an arbitrarily old parent.

The reference implementation TreeTime v0 <a id="cite-2a"></a>[Sagulenko et al. 2018](https://doi.org/10.1093/ve/vex042) [[2](#ref-2)] already constructs slope-based tails (linear in negative-log, exponential in probability) from the outermost trusted points inside its FFT convolution [packages/legacy/treetime/treetime/node_interpolator.py#L231-L256](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/legacy/treetime/treetime/node_interpolator.py#L231-L256). v0 is not the authority for v1 here -- the representations differ (dense FFT vs grid interpolants with a per-side tail policy) -- but its slope-based tail is independent evidence that the flat floor is a simplification, not the intended physics.

### Measured impact

A sensitivity study on `flu/h3n2/20` with `--keep-root --clock-rate=0.003` isolated the effect of the far tail by widening the branch-length grid extent tenfold (`5 * center -> 50 * center`, cap `5.0 -> 50.0` substitutions/site) while holding the grid spacing $\Delta x$ constant (grid points $300 \to 3000$), so the previously flat tail region is filled with the actual computed, decaying likelihood values. A separate control raised the point count to $3000$ at the original extent to isolate resolution from extent.

| Variation                                     | Max internal-node date shift | Internal nodes moved ($>10^{-6}$ y) |
| --------------------------------------------- | ---------------------------- | ----------------------------------- |
| Resolution only (10x points, same extent)     | $3.7 \times 10^{-2}$ y       | 17/17                               |
| Extent (far tail filled with decaying values) | $2.59 \times 10^{1}$ y       | 17/17                               |

The effect is attributable to grid extent, not resolution, by roughly three orders of magnitude. All 36 nodes remained dated in every run; only the 17 internal (non-tip) nodes move, because tip dates are fixed. The deepest node shifted from $1985.7$ to $2011.6$ (NODE_0000016, $25.9$ y), and the direction is consistent across nodes: the flat `Constant` tail biases deep ancestors **older**, and the decaying tail keeps them recent. For H3N2 samples spanning 2000-2013, the decaying result is the physically plausible one.

## Decision

Approved direction: replace the flat far tail with a decaying tail that continues the branch-length likelihood's log-linear slope. Recorded design choices:

- **D1 -- Seam.** Materialize the decay onto appended grid points on the branch-length distribution's long-branch (right) side, at construction, before negation and convolution. The convolution result grid is the sum-support of its operands [packages/treetime-distribution/src/distribution_ops/convolve.rs#L179-L190](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime-distribution/src/distribution_ops/convolve.rs#L179-L190), so extending the source distribution automatically extends the parent message's far-past reach with genuinely convolved decaying values. This is preferred over an evaluation-time `BoundaryBehavior::Exponential` variant on `enum BoundaryBehavior` [packages/treetime-grid/src/grid_fn.rs#L25-L36](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime-grid/src/grid_fn.rs#L25-L36): a bare-grid evaluation policy is never sampled by the convolution and so would not change the parent message's grid.

- **D2 -- Self-terminating extent.** The study showed node dates keep moving as the extent grows, so no fixed extent (and no `Zero` truncation at a fixed extent) is principled. The decay itself supplies the stopping point: append points until the extrapolated density falls below a relative floor $\varepsilon\, f_{\max}$, capped by a physical maximum branch length and the existing $10^{6}$-point safety cap. Below $\varepsilon\, f_{\max}$ the density is numerically negligible and the choice of tail no longer affects results.

- **D3 -- Backward-message tail policy retained.** Keep the backward message's left tail as `Constant` [packages/treetime/src/timetree/inference/backward_pass.rs#L124-L127](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime/src/timetree/inference/backward_pass.rs#L124-L127), now anchored at the extended, negligibly small boundary value rather than at $\mathcal{L}(5 \cdot \text{center})$. This preserves the coverage guarantee that motivated the `Constant` tail (disjoint sibling grids still overlap via the tail, so combining messages does not collapse to `Empty` under `--keep-root`) while removing the bias, because the flat portion now sits at $\le \varepsilon\, f_{\max}$ and contributes negligibly. Switching the tail to `Zero` is the alternative; it is rejected here because it can re-introduce empty products for weakly constrained nodes and would lean on the undated-node warning path [packages/treetime/src/timetree/inference/forward_pass.rs#L52](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime/src/timetree/inference/forward_pass.rs#L52) as a safety net.

### The extrapolation

Given the outermost trusted grid points $t_0 = t_{n-1}$ (boundary) and an interior anchor $t_1 = t_{n-1-m}$, the tail continues the local log-linear slope:

$$f(t) = f(t_0)\,\exp\!\left(\frac{\ln f(t_1) - \ln f(t_0)}{t_1 - t_0}\,(t - t_0)\right), \qquad t > t_0,$$

where $f$ is the plain-probability branch-length density on the grid, $t_0$ is the grid's right boundary, $t_1$ is the interior anchor $m$ points inside the boundary, and $m = \min(3, \lfloor n/3 \rfloor)$ is the slope margin ($n$ is the grid point count). The margin $m > 1$ makes the slope robust to single-point roundoff. The extension is applied only when the estimated slope is negative (a decaying tail); a non-decaying boundary slope leaves the grid unextended.

## Implementation

- **Add a tail-extension operation** on `DistributionFunction` (plain-probability axis), e.g. `fn DistributionFunction.extend_log_linear_tail(side, rel_floor, max_extent)`. It estimates the boundary slope from the outermost $m = \min(3, \lfloor n/3 \rfloor)$ points, and, when the slope is negative, appends grid points at the existing spacing $\Delta x$ continuing the log-linear decay until the extrapolated value drops below `rel_floor * f_max` or the appended axis reaches `max_extent`, whichever comes first. Respect the existing $10^{6}$-point cap. Leave the opposite side untouched.
- **Apply the extension** to the long-branch (right) side inside `fn compute_branch_length_distribution()` [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L23-L60](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L23-L60), after building `normalized_prob` and before constructing the `Distribution::Function`. Extend on the branch-length grid, then map to the time axis, so the decay rate is expressed in the same units the likelihood was computed in. The left (short-branch) side stays at the positive floor; branches cannot be negative, so no left extension is needed.
- **Set constants** as named constants with the derivations above: `rel_floor = 1e-10` (probability floor relative to peak), slope margin `min(3, n/3)`, and reuse the existing `MAX_BRANCH_LENGTH` physical cap as `max_extent`. Do not leave them as bare literals.
- **Keep the backward-message tail** assignment as `Constant` left / `Zero` right [packages/treetime/src/timetree/inference/backward_pass.rs#L124-L127](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime/src/timetree/inference/backward_pass.rs#L124-L127); no change required there, per D3.
- **Check for an existing helper** before writing new array code: search `treetime-utils` array helpers (`packages/treetime-utils/src/array/`) for slope/extension utilities; none is expected to fit, but confirm before adding the method. Stay within `ndarray` operations for the appended values (no drop to `Vec` and back).

## Validation

- **Unit** (`test`): on a constructed exponential density, assert the appended points equal the analytic decay $f(t_0)\exp(s(t-t_0))$ within tolerance; assert termination at `rel_floor * f_max`; assert a non-decaying (flat or increasing) boundary slope leaves the grid unchanged; assert the left side is never extended. Include an asymmetric case so extending the wrong side cannot pass accidentally.
- **Property** (`test-property`): tail values are strictly decreasing; extended density is finite and integrable (finite trapezoidal sum); re-extending an already-decayed grid is approximately idempotent; the analytic decay is matched across random slopes and spacings. Property structural assertions from `treetime-utils` (`prop_assert_array_positive!`, `prop_assert_array_finite!`) apply.
- **Regression / CLI** (smoke): `timetree` on `flu/h3n2/20 --keep-root --clock-rate=0.003` keeps all 36 nodes dated and places the deep internal nodes near the study's extended-grid values (deepest node near $2011$, not $1985$); NODE_0000016 within tolerance of $2011.6$. Commit the assertion referencing the study.
- **Oracle note.** No golden-master test against v0: v0 applies its slope tail at a different seam (FFT convolution result, not the branch-length distribution) and is explicitly not the parity authority here. The analytic exponential decay and the extent self-consistency established by the study are the two independent verification methods.
- **See-red.** Fault-inject each new test (wrong slope sign, impossible tolerance) to confirm it can fail, then restore.

## Knowledge base updates

- **Move** [kb/proposals/exponential-branch-length-tails.md](../proposals/exponential-branch-length-tails.md) out of proposals; its decision is now made and implemented.
- **Update** [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md): the "backward messages have `Constant` left" rationale (the parent is "arbitrarily far in the past ... flat (uninformative)") is superseded for the far region. Record that the branch-length distribution now carries a materialized log-linear decay tail on the long-branch side, that the `Constant` backward-message tail is retained but anchored at a negligible boundary value, and cite the measured node-date impact.
- **Update** [kb/algo/timetree-extrapolation-and-time-clamping.md](../algo/timetree-extrapolation-and-time-clamping.md) cross-references to reflect the decay tail.
- **Delete this ticket** and its source proposal linkage on full resolution, per the ticket lifecycle.

## Related issues

- Source: [kb/proposals/exponential-branch-length-tails.md](../proposals/exponential-branch-length-tails.md)
- Related: [kb/issues/M-timetree-gm-runner-missing-internal-times.md](../issues/M-timetree-gm-runner-missing-internal-times.md) -- internal-node coverage under disjoint multiplication, the failure mode the retained `Constant` tail guards.
- Safety net for the rejected `Zero`-tail alternative in D3: the implemented undated-node warning in `fn propagate_distributions_forward_slot()` [packages/treetime/src/timetree/inference/forward_pass.rs#L52](https://github.com/neherlab/treetime/blob/6ae993ad58d339816df61a243f00aca8b6983041/packages/treetime/src/timetree/inference/forward_pass.rs#L52).

## Glossary

1. <a id="gloss-1"></a> **Backward message.** In the timetree backward pass, the distribution a child sends to its parent, computed as the child's time distribution convolved with the negated branch-length distribution. It answers "when could the parent have existed, given this child?" [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Far tail.** The region beyond a discretized distribution's grid boundary, where the value is supplied by an extrapolation policy (`Error`, `Zero`, `Constant`, or the decay introduced here) rather than by a stored grid point. [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **Log-linear extrapolation.** Continuing a density beyond its grid by extending the straight-line slope of its logarithm, equivalent to exponential decay in probability. [↩](#gloss-use-3)

## References

1. <a id="ref-1"></a> Felsenstein, Joseph. 1981. "Evolutionary trees from DNA sequences: a maximum likelihood approach." _Journal of Molecular Evolution_ 17:368-376. https://doi.org/10.1007/bf01734359 [↩](#cite-1a)
2. <a id="ref-2"></a> Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: maximum-likelihood phylodynamic analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-2a)
