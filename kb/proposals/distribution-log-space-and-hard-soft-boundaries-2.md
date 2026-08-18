# Amendment: convolution grids, backward-pass operations, and mass scope

This document specifies three corrections to the timetree time-distribution design: how convolution against a point or a range behaves, how the backward pass builds its grids, and where mass sizing is used. It amends [distribution-log-space-and-hard-soft-boundaries.md](distribution-log-space-and-hard-soft-boundaries.md), which covers the surrounding design (the switch to log-space storage, the hard/soft boundary framework, and the inference failure that motivated them).

## Background

Terms used throughout:

- **Neg-log storage**: a time distribution stores $-\ln p$ at each grid point. Multiplying distributions is adding ordinates; the likelihood peak is the _least_ ordinate.
- **Grid plus boundary laws**: a distribution is an explicit grid over a finite window, with a boundary law on each side that continues it past the window.
  - **Hard boundary**: probability is zero beyond it; the domain ends there. Near it the density is a power law, $p \propto (t - t_{\text{hard}})^{b}$ (type `HardApproach`, exponent $b$). $b = 0$ is a finite step, so a mode can sit exactly on the edge.
  - **Soft boundary**: an analytic tail continues beyond, log-linear in $t$ (type `SoftTailLaw`, slope $k$: $-\ln p = C + k\,(t - t_{\text{edge}})$).
  - **`Error`**: a side that must never be evaluated.
- **Branch-length factor**: for a branch with $n$ mutations, $L(t) \propto (\mu t)^{n} e^{-\mu t}$, with a hard bound at $t = 0$ and its peak at the maximum-likelihood branch length.
- **Passes**: timetree inference runs a backward pass (leaves to root, building each node's outgoing message) and a forward pass (root to leaves). A coalescent prior adds a per-node weight.
- **Mass windowing** (re-window): trim a distribution's grid to a target fraction of its probability mass, the boundary laws carrying the remainder.

## Convolution against a point or a range

A branch-length distribution convolved with a leaf date behaves differently for an exact date and for an uncertain one:

- **branch ⊛ point** (exact date): the branch distribution shifted along the time axis. Shape and tails do not change; only the position moves.
- **branch ⊛ range** (date interval): the same tails, with the interior washed out. Convolving with a uniform interval smooths the bulk near the peak but leaves each tail's slope intact.

In both cases the tails survive the convolution, so the hard-approach and soft-tail laws carry through directly. These are analytic convolutions: they do not go through the FFT and do not need mass sizing.

## The backward pass has two grid operations

The child fold and the edge convolution build their grids by different rules. They are separate steps.

**Child multiplication** produces the product of the child messages at a node:

1. Multiply all children on one common grid.
2. Take the grid extents from the operands.
3. Take the grid spacing (pitch) from the finest operand.
4. Multiply by evaluating each operand on that grid and adding ordinates (neg-log).
5. Add the coalescent contribution on the same grid.

This step does not use mass sizing. The grid comes from the operand extents and the finest operand's pitch.

**Convolution** propagates a message across an edge:

1. Pick the output grid by mass.
2. Convert to linear (plain) space.
3. Convolve.
4. Convert back to log (neg-log) space.
5. Refit the tails.

Mass sizing selects the output grid here, once per convolution.

## Where mass sizing is used

Mass sizing is not part of branch-length distribution construction. The branch factor has a known closed form, so its grid follows from that form (hard boundary at $t = 0$, peak at the maximum-likelihood branch length, soft tail past the peak). No mass step is needed to build it.

Mass sizing has exactly two consumers, served by the same machinery:

- **Regridding**: choosing the convolution output grid. This is the only regridding in the pass, done once per convolution, not after every multiply.
- **HPD (confidence-interval) reporting**: an interval is a statement about mass, so its integral uses the same mass computation.

This narrows the earlier design, which applied one mass rule at construction and after every operation. Construction is excluded, and the child multiplication sizes its grid from operand extents and the finest pitch rather than from mass.

---

## Extracted knowledge inventory

Each item carries provenance flags:

- 🧮 derived (mathematical derivation)
- 🔬 numerically validated against an analytic oracle
- 💻 traced from current source (verify before relying)
- 📄 from the base proposal
- 🚩 open decision, unverified, or contradicts current code

### Notation

- $t$: time variable (branch length or node time)
- $f, g, h$: probability densities (convolution operands and result)
- $H$: running integral (antiderivative) of $h$
- $L$: branch-length likelihood
- $\mu$: substitution rate (expected mutations per unit time)
- $n$: mutation count on a branch
- $p, q$: power-law orders of two operands near a shared edge
- $w$: half-width of a date range
- $s$: center of a date range
- $A$: amplitude of a point or range operand
- $a$: location of a hard boundary
- $b$: hard-approach exponent (power-law order at a hard boundary)
- $k$: soft-tail log-slope
- $C$: a constant factor
- $t_0$: position of a point operand (a fixed date)
- $t_{\text{hard}}$: hard-boundary location of an operand
- $x_{\min}, x_{\max}$: grid lower and upper extent
- $dx$, $dx_1$, $dx_2$: grid spacing (pitch); the two operand spacings
- $N_1, N_2$: the two operand grid point counts
- $\epsilon$: tail mass fraction left outside the grid, per side
- $\ast$ (written ⊛ in prose): convolution
- $\delta$: Dirac delta, a unit point mass; $\mathbf{1}_{[\alpha,\beta]}$: indicator (uniform box) on $[\alpha,\beta]$
- $\alpha, \beta$: generic interval endpoints; $L_1, L_2$: widths of two boxes

### Convolution mathematics

Convolution combines a branch-length distribution with a leaf-date operand. In continuous form,

$$(f \ast g)(t) = \int_{-\infty}^{\infty} f(\tau)\,g(t - \tau)\,d\tau,$$

evaluated in plain (linear) probability space, then stored back in neg-log. The three operand kinds are a **point** (an exact date, a Dirac delta $A\,\delta(t - t_0)$), a **range** (an uncertain date, a uniform box), and a **function** (a full gridded distribution, such as a branch-length curve).

Two facts hold for every case:

- **Support is the Minkowski sum.** If $f$ is supported on $[\alpha_1, \beta_1]$ and $g$ on $[\alpha_2, \beta_2]$, then $f \ast g$ is supported on $[\alpha_1 + \alpha_2,\ \beta_1 + \beta_2]$; a convolution of compactly supported functions is compactly supported.
- **Each edge transforms on its own.** A hard edge (density reaches zero) and a soft edge (analytic tail) are reshaped independently by the other operand. Under convolution a one-sided power law $(t - t_{\text{hard}})^{b}$ gains one order, becoming $(t - t'_{\text{hard}})^{b+1}$, because the kernel integrates it once (Cauchy's repeated-integration formula). Under multiplication the exponents add with no gain.

Each case below is a self-contained spec: operands, output support, bulk (interior) computation, boundary laws per side, and implementation steps.

<details>
<summary>C1. point ⊛ point 🧮</summary>

**Operands.** Two point masses, $A_1\,\delta(t - t_1)$ and $A_2\,\delta(t - t_2)$.

**Result and support.** The delta is the convolution identity ($f \ast \delta_{t_0} = f(\cdot - t_0)$), so

$$\big(A_1\,\delta_{t_1}\big) \ast \big(A_2\,\delta_{t_2}\big) = A_1 A_2\;\delta\!\big(t - (t_1 + t_2)\big).$$

The result is a single point at $t_1 + t_2$: positions add, amplitudes multiply.

**Bulk.** None; a point has no interior. In neg-log the ordinate adds: $-\ln(A_1 A_2) = -\ln A_1 - \ln A_2$.

**Boundary laws.** A point is a zero-width hard support; both sides are hard at $t_1 + t_2$. No tails.

**Implementation.** Emit a point $(t_1 + t_2,\ A_1 A_2)$; no grid is allocated.

</details>

<details>
<summary>C2. point ⊛ range 🧮</summary>

**Operands.** A point $A_p\,\delta(t - t_0)$ and a range, the uniform box $A_r\,\mathbf{1}_{[r_1, r_2]}$.

**Result and support.** The delta translates the box:

$$\big(A_p\,\delta_{t_0}\big) \ast \big(A_r\,\mathbf{1}_{[r_1, r_2]}\big) = A_p A_r\;\mathbf{1}_{[r_1 + t_0,\ r_2 + t_0]}.$$

Support $[r_1 + t_0,\ r_2 + t_0]$: the range shifted by $t_0$, width unchanged.

**Bulk.** A uniform box of constant height $A_p A_r$; the neg-log ordinate is the constant $-\ln(A_p A_r)$.

**Boundary laws.** Both edges are hard (the box steps to zero at each end); no soft tail.

**Implementation.** Shift both range endpoints by $t_0$, scale the height by $A_p$, declare both sides `Hard`.

</details>

<details>
<summary>C3. range ⊛ range 🧮</summary>

**Operands.** Two boxes $A_1\,\mathbf{1}_{[\alpha_1, \beta_1]}$ and $A_2\,\mathbf{1}_{[\alpha_2, \beta_2]}$, of widths $L_1 = \beta_1 - \alpha_1$ and $L_2 = \beta_2 - \alpha_2$.

**Result and support.** The convolution of two boxes is a trapezoid, degenerating to a triangle when $L_1 = L_2$. Its value at $t$ is the length of the overlap of $[\alpha_1, \beta_1]$ with the reflected, shifted second box $[t - \beta_2,\ t - \alpha_2]$:

$$(f \ast g)(t) = A_1 A_2 \cdot \big|\,[\alpha_1, \beta_1] \cap [t - \beta_2,\ t - \alpha_2]\,\big|.$$

Support $[\alpha_1 + \alpha_2,\ \beta_1 + \beta_2]$, width $L_1 + L_2$; peak height $A_1 A_2 \min(L_1, L_2)$.

**Bulk.** Piecewise-linear, with breakpoints at $\alpha_1 + \alpha_2$, then $+\min(L_1, L_2)$, then $\beta_1 + \beta_2 - \min(L_1, L_2)$, then $\beta_1 + \beta_2$: it rises with slope $A_1 A_2$, plateaus (a triangle has no plateau), then falls with slope $-A_1 A_2$.

**Boundary laws.** Both outer edges are hard (compact support). The two inner breakpoints are slope discontinuities (corners), not boundaries. No soft tail.

**Implementation.** Emit piecewise-linear nodes at the four breakpoints; declare both ends `Hard`.

</details>

<details>
<summary>C4. point ⊛ function 🧮</summary>

**Operands.** A point $A\,\delta(t - t_0)$ and a gridded density $f$ with its boundary laws.

**Result and support.** Pure translation and scaling:

$$g(t) = A\,f(t - t_0).$$

Support is $f$'s support shifted by $t_0$; width and shape unchanged.

**Bulk.** Same grid spacing $dx$ and count $N$, with $x_{\min} \to x_{\min} + t_0$. In neg-log the ordinates shift uniformly: $y_g(t) = y_f(t - t_0) - \ln A$.

**Boundary laws** (both carried, no refit):

- **Hard side.** `HardApproach` with exponent $b$ unchanged; $t_{\text{hard}} \to t_{\text{hard}} + t_0$.
- **Soft side.** `SoftTailLaw` with slope $k$ unchanged. The law is edge-relative (stored against the edge ordinate), and the edge ordinate shifts by the same $-\ln A$, so it stays consistent without refitting.

**Implementation.** Shift the grid origin by $t_0$, add $-\ln A$ to all ordinates, shift $t_{\text{hard}}$, copy both laws.

</details>

<details>
<summary>C5. range ⊛ function 🧮🔬</summary>

**Operands.** A range (box of half-width $w$, center $s$, height $A_r$) and a gridded density $f$. Factor the range as a point at $s$ (handled by C4) followed by a zero-centered box of half-width $w$. Let $h$ be $f$ after the C4 shift-and-scale; the box convolution acts on $h$.

**Bulk (cumulative-integral difference).** Convolving with the box is a moving average over a window of width $2w$; it blurs the interior. Compute it as a difference of running integrals, not a masked point sum:

$$g(t) = \int_{t - w}^{t + w} h(u)\,du = H(t + w) - H(t - w), \qquad H(x) = \int_{-\infty}^{x} h(u)\,du.$$

Build $H$ **hard-anchored**:

- $H(x) = 0$ for $x \le t_{\text{hard}}$, with no grid nodes below the bound. Integrating a zero node below the bound against the finite boundary node fabricates a triangle of mass that does not exist.
- Over the support, $H$ is the trapezoidal running integral of the piecewise-linear $h$.
- Beyond a soft edge, extend with the closed-form tail integral $H(x) = H(x_{\max}) + \int_{x_{\max}}^{x} h$ using the log-linear law.

**Convergence.** Because $h$ is piecewise-linear, $H$ is piecewise-quadratic, and linear interpolation of $H$ has $O(dx^2)$ error. The rectangle-mask sum $dx \sum_{|t_j - t| \le w} h_j$ mis-weights the partial end cell at the hard turn-on by $O(1)$ and does not converge under $dx$ refinement.

**Support.** Widened by $w$ on each side: $[\,t_{\text{hard}} - w,\ x_{\max} + w\,]$, then translated by the center $s$.

**Boundary laws.**

- **Hard side.** The box integral of $(u - t_{\text{hard}})^{b}$ over the partial window is $\propto (t - (t_{\text{hard}} - w))^{b+1} / (b + 1)$: the exponent gains one and the bound moves outward by $w$. A step ($b = 0$) becomes a linear ramp ($b = 1$) spanning $[t_{\text{hard}} - w,\ t_{\text{hard}} + w]$; $b = 1$ becomes quadratic. Declare `HardApproach` with exponent $b + 1$ at $t_{\text{hard}} - w$.
- **Soft side.** The window integral of $e^{-k u}$ is $\int_{t-w}^{t+w} e^{-k u}\,du = e^{-k t}\,\dfrac{2\sinh(k w)}{k}$: the slope $k$ is unchanged, only a constant factor differs. Refit `Linear` (slope $k$) on the edge moved out by $w$.

**Oracle (acceptance test).** For $h(t) = e^{-\mu t}$ on $t \ge 0$ (a hard step at $0$):

$$g(t) = \frac{1}{\mu}\Big[\,e^{-\mu \max(0,\ t - w)} - e^{-\mu (t + w)}\,\Big].$$

Support $[-w, \infty)$; the step becomes a linear ramp on $[-w, w]$; the right-tail slope stays $-\mu$.

**Implementation.** Build hard-anchored $H$; widen the grid by $\lceil w / dx \rceil$ cells per side; evaluate $g(t) = H(t + w) - H(t - w)$; set the left law `HardApproach`($b + 1$) at $t_{\text{hard}} - w$ and the right law `Linear` (slope $k$); translate the result by the center $s$.

</details>

<details>
<summary>C6. function ⊛ function 📄💻🧮</summary>

**Operands.** Two gridded densities $f$ (spacing $dx_1$, $N_1$ points) and $g$ (spacing $dx_2$, $N_2$ points).

**Method (FFT, in plain space).** A convolution is an integral, not a pointwise operation, so it is carried out in plain probability space:

1. **Resample** both to a common spacing $dx = \min(dx_1, dx_2)$; the FFT needs a shared grid and the finer spacing is the resolution floor.
2. **Peak-normalize.** Convert neg-log to plain by subtracting the peak ordinate and exponentiating (prevents underflow); record the peak offset in neg-log units.
3. **Full-mode FFT convolution.** The output has $N_1 + N_2 - 1$ points and already spans the Minkowski sum.
4. **Convert back** to neg-log and restore the peak offset.
5. **Mass-size** the output grid (the FFT grid is large; trim to the target mass fraction).

**Support.** $[\,x^{f}_{\min} + x^{g}_{\min},\ x^{f}_{\max} + x^{g}_{\max}\,]$.

**Bulk.** Trustworthy only above the FFT roundoff floor (about $10^{-16}$ of peak in plain space); below that the output is numerical noise.

**Boundary laws** (reconstructed by fit, never read from the FFT):

- **Hard side.** Exponents add and gain one ($t^{p} \ast t^{q} \sim t^{p + q + 1}$); the result edge is the sum of the operand hard edges. Refit `HardApproach` from the outermost trusted points.
- **Soft side.** Refit `Linear` slope from the outermost trusted points; discard everything below the roundoff floor.

**Implementation.** Resample, peak-normalize, full FFT, denormalize, refit both tails from the trusted band, then mass-window. This is the only case that mass-sizes its own output.

</details>

### Code divergences

Verify these against current source before relying on line numbers.

<details>
<summary>K11. Analytic convolutions drop operand tails 💻</summary>

`convolution_point_function` and `convolution_range_function` build via `from_start_dx_values`, then `from_grid_array`, which **resets both boundary sides to `Error`** [packages/treetime-grid/src/grid_fn.rs#L62-L76](packages/treetime-grid/src/grid_fn.rs#L62-L76). The operand tails are lost.

The backward and forward passes then fit fresh tails **after** each convolution [packages/treetime/src/timetree/inference/backward_pass.rs#L143-L146](packages/treetime/src/timetree/inference/backward_pass.rs#L143-L146), [packages/treetime/src/timetree/inference/forward_pass.rs#L235-L238](packages/treetime/src/timetree/inference/forward_pass.rs#L235-L238). That recovers a soft slope, but it cannot recover a dropped hard bound or truncated support.

</details>

<details>
<summary>K12. convolution_range_function has three defects 💻🔬</summary>

Main path [packages/treetime-distribution/src/distribution_ops/convolve.rs#L130-L161](packages/treetime-distribution/src/distribution_ops/convolve.rs#L130-L161), mask at [packages/treetime-distribution/src/distribution_ops/convolve.rs#L150-L159](packages/treetime-distribution/src/distribution_ops/convolve.rs#L150-L159):

1. **Truncated support**: it reuses the input extent, so there is no widening by $w$.
2. **Dropped tails**: via `from_start_dx_values`.
3. **Biased edges**: the mask sum ignores off-grid mass, so the outermost retained ordinates read low. The later tail-refit then reads these corrupted points.

</details>

<details>
<summary>K13. Branch factor is mass-sized at construction, against the intended design 💻🚩</summary>

`compute_branch_length_distribution` currently **mass-sizes the freshly built factor** [packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L66-L72](packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L66-L72). The intended design builds it analytically from its closed form, with no re-window at construction.

This contradicts landed code and needs an explicit decision.

</details>

<details>
<summary>K14. Branch-length closed form and its analytic grid 📄💻</summary>

The factor carries a **Poisson indel term**:

$$L(t) \propto (\mu t)^{n}\, e^{-\mu t} \times (\text{Poisson indel term})$$

Its grid follows directly from the form:

- hard bound at $t = 0$,
- peak at the ML length,
- soft log-linear tail past the peak.

</details>

<details>
<summary>K15. Mass machinery location 💻</summary>

The mass functions live in [packages/treetime-distribution/src/distribution_ops/mass_domain.rs](packages/treetime-distribution/src/distribution_ops/mass_domain.rs): `total_mass`, `mass_bounded_domain`, `rewindow_to_mass`. `rewindow_to_mass` trims to $\ge 1 - 2\epsilon$ of the mass.

The **same functions** serve regridding and should serve HPD.

</details>

### Confidence intervals (HPD)

<details>
<summary>K16. Confidence integrals ignore tail mass 💻</summary>

`quantile` and `hpd_region` integrate a **grid-only trapezoidal CDF** and ignore the boundary-law tail masses [packages/treetime-distribution/src/distribution_core/distribution.rs#L266-L323](packages/treetime-distribution/src/distribution_core/distribution.rs#L266-L323), [packages/treetime-distribution/src/distribution_core/distribution.rs#L550-L635](packages/treetime-distribution/src/distribution_core/distribution.rs#L550-L635).

Re-window deliberately leaves $2\epsilon$ mass in the tails, so a grid-only integral normalizes against $1 - 2\epsilon$ of the true mass and **biases every reported interval inward**.

**Fix**: integrate the same peak-relative mass profile that `rewindow_to_mass` uses (the grid trapezoid plus both closed-form tail masses).

</details>

<details>
<summary>K17. Related open tickets 💻</summary>

- [kb/issues/M-timetree-confidence-marginal-hpd-disabled-under-neglog.md](../issues/M-timetree-confidence-marginal-hpd-disabled-under-neglog.md)
- [kb/issues/M-timetree-confidence-interval-deficiencies.md](../issues/M-timetree-confidence-interval-deficiencies.md)

</details>

### Backward pass and coalescent

<details>
<summary>K18. Coalescent factor is message-only 📄💻</summary>

The coalescent factor belongs to the **outgoing (backward) message only**, never the stored node distribution. It is added on the node grid before the single end-of-node re-window.

</details>

<details>
<summary>K19. Common grid for the child product 📄</summary>

Child multiplication builds the common grid by **intersecting hard sides to the tightest bound and widening soft sides to the loosest**, with spacing taken from the finest operand.

</details>

### Open decisions and caveats

<details>
<summary>K20. Convolution grid handling is undecided 🚩</summary>

Grid selection and manipulation before and after convolution is not settled. A concrete scheme is proposed and validated but not implemented; per the project scientific-change gate it needs explicit consent.

</details>

<details>
<summary>K21. Which left-law to declare for range ⊛ function 🧮🚩</summary>

- **Default**: `HardApproach` with exponent $b + 1$ at $a - w$ (sub-grid accurate).
- **Plain `Hard` at $a - w$**: acceptable only when $w \gg dx$, so the ramp is grid-resolved.
- **When $w \lesssim dx$**: the ramp is sub-cell, and only the `HardApproach` law represents it.

</details>

<details>
<summary>K22. Grid-map sizing rule 🧮</summary>

- The FFT path self-sizes by mass.
- Every analytic path produces a naturally compact, tail-complete result and defers sizing to the pass re-window.
- The resolution floor is the finest operand $dx$.

</details>

<details>
<summary>K23. branch ⊛ point result is derived, not confirmed 🚩</summary>

The domain experts did not specify the branch ⊛ point case. The pure-shift result in C4 is derived from the convolution, not a confirmed decision. Confirm whether it is the intended design.

</details>
