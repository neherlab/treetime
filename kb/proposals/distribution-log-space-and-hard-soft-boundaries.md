# Log-space time distributions with hard/soft support boundaries

Move timetree time distributions to a logarithmic y-axis and replace the current single-meaning support boundary with an explicit distinction between **hard** boundaries (a fact about the distribution: probability is zero beyond) and **soft** boundaries (a choice of representation: an analytic tail continues beyond). The explicit grid becomes a _sensible domain_ holding a target fraction of the probability mass, re-derived once per operation, rather than the distribution's support.

**Type**: New v1 design. Supersedes part of [kb/decisions/timetree-inference-pass-boundary-tails.md](../decisions/timetree-inference-pass-boundary-tails.md) and amends [kb/decisions/distribution-intersection-grid-resolution.md](../decisions/distribution-intersection-grid-resolution.md).

## Motivating defect

`treetime timetree` on a 967-leaf Ebola-scale dataset (197,209 sites, `--clock-rate 5.7e-5`) reports `log_lh_pos=NaN` and emits an Auspice tree in which no internal node carries `num_date`. The command exits 0 with no warning.

Chain of causation:

1. `create_simple_grid()` builds the branch-length grid over `[0.01 · one_mutation, max(5 · branch_length, 10 · one_mutation)]`. For the near-zero branches that dominate this tree the extent is `10 · one_mutation`. Because the extent is a multiple of `one_mutation`, its width in time is exactly that multiple of mutations: with `μ = clock_rate · L = 11.24` mutations/yr, the distribution spans **0.8887 years**.
2. A leaf's `msg_to_parent = point(date) ⊛ (−branch_dist)` therefore has support `[date − 0.8896, date − 0.0009]`.
3. Samples span 2018–2025. Any cherry or polytomy whose children are sampled more than ~0.89 yr apart produces messages with disjoint supports. Observed at node 603: accumulator `[2021.5501, 2022.4197]`, next child's message `[2024.6816, 2025.5703]`.
4. `multiply_function_function()` intersects raw grid bounds and returns `Distribution::empty()` on `Disjoint`. The `left_extrap = Constant` tail that the backward pass deliberately attached is ignored — multiplication never consults operand tails.
5. `Empty` is stored on the node, its outgoing message is `Empty`, and every ancestor to the root follows. Measured: 4 initial collapses → 58 internal nodes `Empty` including the root → all 471 internal nodes left without a time in each of the 2 iterations.
6. `compute_positional_log_lh()` skips every edge because both endpoint times are `None`, so `count == 0`, the function returns `None`, and the metric prints as `NaN`.

Confirmed by experiment: flooring `max_bl` at 1e-3 subs/site eliminates all 943 collapses, yields `log_lh_pos = -2089.80`, and populates `num_date` on all 1438 nodes (root 2016.91).

## The grid was never the problem

The mass actually contained in the current window, integrating the Gamma likelihood `L(t) ∝ (μt)^n e^{-μt}`:

| branch                              | grid extent | mass outside |
| ----------------------------------- | ----------- | ------------ |
| n=0 (the `10 · one_mutation` floor) | 10 mut      | 4.5e-5       |
| n=1 (`5 · center`)                  | 5 mut       | 4.0e-2       |
| n=2                                 | 10 mut      | 2.8e-3       |
| n=3                                 | 15 mut      | 2.1e-4       |
| n≥5                                 | ≥25 mut     | <1.4e-6      |

For the zero-length branches that dominate this dataset the existing window already holds 99.9955% of the mass. It is a good _sensible domain_; it was only ever wrong as a _support_. Two consequences:

- The fix is boundary semantics, not a wider grid. [kb/tickets/timetree-switch-branch-grid-to-nonuniform-spacing.md](../tickets/timetree-switch-branch-grid-to-nonuniform-spacing.md) drops from correctness to accuracy.
- The `5 ·` peak multiple is the one genuinely tight case (4% of mass outside at n=1). It should be replaced by the same mass criterion used for re-windowing, not tuned.

## v0 reference

v0 already does both halves of this proposal, though not under these names.

- **Log axis**: v0 stores neg-log probability throughout (`packages/legacy/treetime/treetime/distribution.py:216`).
- **Analytic tails**: `NodeInterpolator.convolve_fft()` (`packages/legacy/treetime/treetime/node_interpolator.py:231-256`) constructs slope-based tails from the outermost trusted points — linear in neg-log, i.e. exponential in probability — extending them only where the tail decays away from support.
- **Wide support**: `BranchLenInterpolator` (`packages/legacy/treetime/treetime/branch_len_interpolator.py:36-60`) always reaches `MAX_BRANCH_LENGTH`, so v0 supports are never disjoint in practice.

v1 kept the finite grid but dropped both the log axis and the analytic tails. This proposal restores them; it reduces divergence from v0 rather than adding to it.

## Core reframing

> A **hard** boundary is a fact about the distribution. It is immovable and probability is zero beyond it.
> A **soft** boundary is a choice of representation. It is freely movable and an analytic law continues beyond it.

Everything else follows mechanically. Multiplication intersects the facts and re-chooses the representation. Re-windowing never _drops_ a tail, it re-derives one, so no mass leaks however many operations a distribution passes through.

Critically, **a hard boundary is not necessarily a zero of the density**. For a zero-mutation branch `L(t) = e^{-μt}` is _maximal_ at `t = 0`. In the node-603 collapse the correct answer places the parent exactly at its earlier child's date — the mode sits on the hard boundary. Any rule that assumes "hard ⇒ density → 0", or that trims a hard edge inward, destroys precisely the case this proposal exists to fix.

## Part A — logarithmic y-axis

Switch timetree time distributions from `Plain` to `NegLog` (`packages/treetime-distribution/src/policy.rs`). Both policies already exist and `Distribution<Y>` is already generic over them.

### Why

1. **The tail stops being a special case.** In log space the Gamma tail `n·log t − μt` is smooth and asymptotically straight, so the existing piecewise-linear grid is already the right basis for it. "Exponential extrapolation" reduces to "keep interpolating linearly past the edge" — a `BoundaryBehavior::Linear` variant on `GridFn`, not new distribution machinery.
2. **Multiplication is addition.** Exact, associative, no intermediate normalization.
3. **Dynamic range.** The value that decides node 603 is 9.1e-12 of peak; the plateau of a true GTR likelihood is ~`e^{-273000}` below peak. Neither is a problem in log space.
4. **Boundary laws compose in closed form** (see Part B).

### What it closes

- Backward-pass child fold: the per-step `normalize()` existed only to stop `(peak)^N` plain-space underflow. Under `NegLog` the fold is addition, so the underflow and the per-step normalize both disappear.
- [M-distribution-plain-division-fixed-floor.md](../issues/M-distribution-plain-division-fixed-floor.md) — division is subtraction; `NegLog::safe_divisor` is already the identity.
- [N-distribution-function-product-negative-roundoff.md](../issues/N-distribution-function-product-negative-roundoff.md) — linear interpolation of a log-density cannot produce a negative amplitude, so the `~-1e-26 → ln → NaN` path cannot recur.

The coalescent contribution also stops round-tripping: `distribution_multiply_by_fn` already computes in neg-log and becomes plain addition.

### Prerequisite

`Distribution::likely_time()` must select the likelihood peak per policy: under `NegLog` the peak is the _least_ stored ordinate, not the greatest. Every inference node depends on this once the switch lands. This prerequisite landed: `likely_time` is now policy-aware (`function.rs`).

### Normalization

Under `NegLog`, normalization is _subtraction of the minimum ordinate_ (peak → 0): a shift, exactly representable, no scaling. `max_value()` / `scale_by()` / `normalize()` (`distribution.rs:214-243`) need policy-aware replacements.

### Convolution stays in plain space

`SupportsConvolution` is implemented only for `Plain`, correctly — an integral is not a pointwise operation. Each convolution becomes log → peak-subtracted exp → FFT → log, once per edge per pass. See Part D for the accuracy constraint this imposes.

## Part B — hard/soft boundary semantics

### Classification

| Behavior                | Class      | Meaning                                                                                          |
| ----------------------- | ---------- | ------------------------------------------------------------------------------------------------ |
| `Zero` (rename: `Hard`) | hard       | Domain terminates. p = 0 beyond. Restricts the result domain.                                    |
| `Linear` (new)          | soft       | Log-linear law continues beyond. Extends the result domain.                                      |
| `Constant`              | soft       | Flat tail. Retained only for genuinely uninformative edges.                                      |
| `Error`                 | undeclared | Restricts, and evaluation beyond is a programming error. Must never be silently treated as soft. |

`Constant` should be retired from the inference passes. It is non-integrable, so it silently corrupts `quantile()` and `hpd_region()` — the confidence-interval path. `Linear` is integrable in closed form.

### Multiplication rule

For each side independently:

```text
hard side:  result bound = tightest (innermost) hard bound among operands
soft side:  result bound = loosest  (outermost) soft bound among operands
mixed:      a hard bound always dominates a soft bound on the same side
```

Each operand is then evaluated over the result domain, using its declared law wherever the query falls outside its own grid. Division now uses this same per-side rule (`fn divide_function_by_function()`; see [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md)): it is no longer asymmetric with a privileged dividend, and differs from multiplication only in combining ordinates by subtraction and refitting its soft result tails from the combined grid.

Applied to node 603:

| side  | accumulator `[2021.5501, 2022.4197]` | message `[2024.6816, 2025.5703]` | result                   |
| ----- | ------------------------------------ | -------------------------------- | ------------------------ |
| right | hard (child's date)                  | hard                             | tightest → **2022.4197** |
| left  | soft                                 | soft                             | loosest → **2021.5501**  |

Non-empty. Both operands decay leftward, so the product peaks at the hard right edge: the parent is placed at its earlier child's date and the 2025 child absorbs a ~3.15 yr branch. That is the correct MAP for data that genuinely conflicts with a 5.7e-5 clock.

### Boundary laws compose exactly under multiplication

Under `NegLog` both laws are linear in a transformed coordinate:

- soft tail: `y = a + k·(t − t_edge)`, linear in `t`
- hard approach: `y = a − b·log|t − t_hard|`, linear in `log|t − t_hard|` (Part C)

Multiplication is addition, and the sum of two linear functions is linear in the same coordinate. **Multiplication therefore propagates both laws in closed form, with no refitting.** Only convolution requires a refit.

### Invariant

`Distribution::Empty` becomes reachable **only** from genuinely disjoint hard domains. It must never arise from numerical collapse. Given that silent `Empty`-poisoning produced the original defect, this should be a checked invariant, not a consequence.

### Policy must survive regridding

`GridFn::from_grid_array()` resets both sides to `Error` (`grid_fn.rs:67-80`), so every path through `scale_y` / `normalize` currently discards the tails. In the backward child fold the accumulator is normalized after every multiply, so the new rule would never fire from the second child onward. Regridding must carry the policy _and_ the fitted law. `resample()` already does this (`grid_fn.rs:485-497`); `from_grid_array` callers do not.

## Part C — interpolation between the first grid point and a hard boundary

**This region needs its own law, distinct from the soft-tail law, and it is the subtlest part of the design.**

### The problem

Near a hard boundary the density behaves as a power law: for a branch with `n` mutations, `p(t) ∝ t^n` as `t → 0⁺`. In neg-log:

```text
y(t) = −n·log t + μt + C
```

For `n ≥ 1` this **diverges logarithmically** as `t → t_hard`. Three failure modes follow from ignoring it:

1. **Storing the boundary ordinate.** `y(t_hard) = +∞` for `n ≥ 1`. Any linear interpolation into a cell with an infinite endpoint returns `+∞` across the whole cell, wrongly zeroing a finite region of density. `to_plain(+∞) = 0` hides this rather than erroring.
2. **Linear-in-`t` interpolation across the first cell.** The true `y` is convex-divergent, so linear-in-`t` is a poor approximation exactly where the density is changing fastest.
3. **Assuming the ordinate diverges at all.** For `n = 0`, `y(t) = μt + C` is _finite_ at the boundary and the density is _maximal_ there. The edge is a genuine step discontinuity, not a zero. This is the common case for the zero-length branches in the motivating dataset, and it is the case whose mode must be preserved exactly.

### Proposed law

Fit a two-parameter approach law from the innermost grid points, by least squares on `(log(t − t_hard), y)`:

```text
y(t) = a − b·log(t − t_hard),   b ≥ 0
```

- `b = 0` recovers the finite/step case (`n = 0`) exactly.
- `b > 0` recovers the power-law vanishing case (`n ≥ 1`) exactly.

One law covers both; no branching on `n`, which is not available after convolution anyway. Constrain `b ≥ 0`: a negative fit means the density _increases_ into the boundary faster than any power law, which is unphysical here and indicates noise in the innermost points — clamp to `b = 0`.

### Grid placement rule

> Never store an infinite ordinate on a grid point.

- When `b > 0` (divergent), the grid stays strictly inside the hard boundary and the approach law covers `[t_hard, t_first)`. `create_simple_grid`'s existing `min_bl = 0.01 · one_mutation` floor already does this for the `t = 0` edge; formalize it as a rule rather than an incidental guard for the Poisson indel term.
- When `b = 0` (finite), the hard edge is stored as an exact grid endpoint, preserving the step discontinuity and any mode sitting on it. This is compatible with the existing spacing contract, which already includes both analytical endpoints.

### Mass in the approach region

Closed form, mirroring the soft tail's `exp(−a)/k`:

```text
∫_{t_hard}^{t_1} exp(−y) dt = exp(−a) · (t_1 − t_hard)^{b+1} / (b + 1)
```

Both boundary regions therefore contribute analytically to the mass integral used by Part D, and neither needs to be discretized.

### Composition

- **Multiplication**: exponents add (`b_result = b_a + b_b`), which is exactly the addition of the two log-space laws. Closed form, no refit — consistent with Part B.
- **Convolution**: exponents add _and gain one_ (`t^m ⊛ t^n ∼ t^{m+n+1}` near the combined edge), and the edge itself is the sum of the operand edges. This must be refit from the result, not propagated.
- **Quantile / HPD**: integration must not smooth across a `b = 0` step. The discontinuity is real.

## Part D — adaptive sensible domain

The explicit grid represents a target mass fraction (>0.999 proposed); the boundary laws carry the rest. After each operation the domain is re-derived and the grid resampled to fit.

```text
domain = [ max(hard_lo, q(ε)),  min(hard_hi, q(1−ε)) ]
```

with the ε-trim **suppressed on any side where the mode sits on the hard bound** — otherwise node 603's answer is trimmed away.

Mass-based (quantile) rather than magnitude-based (`y < δ · y_max`): mass is the meaningful invariant, and both boundary regions contribute in closed form, so the integral is cheap. The same criterion replaces the ad-hoc `max(5 · center, 10 · one_mutation)` extent in `create_simple_grid()` — one rule, applied at construction and after every operation.

### Backward pass restructuring

Pointwise addition requires co-located grids, which changes the shape of the child fold:

1. Choose the node's working grid once — hard bounds intersected, soft extent chosen generously, spacing from the finest operand.
2. Resample each child message onto it exactly once.
3. Sum, in any order.
4. Re-window and resample once, at the end.

Addition is exact and associative with no intermediate normalization, so each message is resampled exactly once regardless of fan-out. Node 880 (25 children) goes from 24 intersect + resample + normalize rounds to 25 resamples and a sum. This removes the resampling-drift concern entirely and makes a balanced multiplication tree unnecessary.

### Point budget

Fixed `n` with adaptive `dx`, plus a resolution floor so that multiplying a narrow distribution by a wide one cannot ratchet `dx` coarser than the narrow operand had. The `MAX_GRID_POINTS = 1_000_000` ceiling stops being load-bearing once domains are mass-bounded rather than support-bounded.

## Design axes

| #   | Axis                  | Options                                                                               | Recommendation                                                                  |
| --- | --------------------- | ------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| 1   | Domain criterion      | mass-based ε-quantile / magnitude-based `δ·y_max`                                     | mass-based; tail mass is closed-form on both sides                              |
| 2   | ε                     | 1e-3 / 1e-4 / configurable                                                            | single named constant, 5e-4 per side                                            |
| 3   | Tail slope estimation | last two points / least squares over outermost window / analytic seed at construction | least-squares window, with decay guard; analytic seed where the family is known |
| 4   | `Constant` retention  | retire entirely / keep for non-inference uses                                         | retire from inference passes; non-integrable                                    |
| 5   | `Zero` naming         | keep / rename to `Hard`                                                               | rename — under this design it is a domain declaration, not a written y-value    |
| 6   | Point budget          | fixed n / adaptive n                                                                  | fixed n with a `dx` floor                                                       |
| 7   | Empty invariant       | documented / debug-asserted / checked error                                           | checked error; this is the original defect                                      |

## Coupling

`YAxisPolicy::supports_zero_boundary()` rejects a `Zero` tail under `NegLog`, because literal `0.0` is the multiplicative identity there. That guard alone blocks Part A. Under Part B "hard" is a _domain_ declaration and no value is written outside the grid, so the guard has nothing left to protect and is removed. **Parts A and B must land together; neither works alone.**

## Accuracy constraint on convolution

The FFT operates in plain space, so its output is trustworthy only in the bulk — anything below ~1e-16 of peak is roundoff. The convolved message's soft tail and hard approach law must both be **reconstructed by fit in log space from the outermost trusted points**, not read out of the FFT. This is why v0 builds tails inside `convolve_fft` rather than carrying a policy through it, and it is the single most important detail to get right: node 603's decisive value is 9.1e-12 of peak, well inside the region the FFT cannot produce but a log-linear fit reproduces exactly.

Reproduce v0's guard: extend a tail only where it decays away from support; otherwise clamp.

Note also that log-linear extrapolation from the grid edge slightly _over_-estimates the far tail, because the true log-slope `n/t − μ` steepens toward `−μ`. The error is conservative — it can never manufacture a spurious zero — and can be removed where the family is known by seeding the law analytically at construction.

## Sequencing

1. Policy-aware `likely_time` -- prerequisite, landed.
2. `BoundaryBehavior::Linear` + hard/soft classification on `GridFn`; policy and fitted law survive `from_grid_array` / `normalize`.
3. Hard-boundary approach law (Part C) with the grid-placement rule.
4. Multiplication rule (Part B) + common-grid addition in the backward pass.
5. Switch timetree distributions to `NegLog`; remove `supports_zero_boundary`.
6. Convolution round-trip with post-hoc law refitting.
7. Adaptive sensible domain (Part D), replacing `create_simple_grid`'s ad-hoc extent.
8. `Empty`-only-from-disjoint-hard invariant.

Steps 2–5 are one atomic behavioral change; 1 precedes them, 6–8 follow.

## Validation plan

- **Regression**: the motivating dataset must produce a finite `log_lh_pos` and `num_date` on all 1438 nodes. Currently `NaN` and 967/1438.
- **Golden masters**: `test_gm_runner_marginal_dense`, `test_gm_runner_marginal_sparse`, `test_gm_runner_poisson` are `#[ignore]`d at `epsilon = 1e-6`. Un-ignoring them is the acceptance criterion already recorded in [M-timetree-branch-grid-uniform-resolution.md](../issues/M-timetree-branch-grid-uniform-resolution.md).
- **Unit — disjoint supports**: two messages with disjoint grids multiply to a non-`Empty` result whose mode is on the tighter hard bound.
- **Unit — tail accuracy**: fitted log-slope matches the analytic Gamma tail to a stated tolerance 2–3 grid widths out.
- **Unit — hard approach**: `b = 0` preserves a mode on a hard edge exactly; `b > 0` matches `t^n` to tolerance across the first cell; no `+∞` is ever stored.
- **Unit — mass conservation**: repeated re-windowing of a fixed distribution conserves total mass and does not drift the mode.
- **Unit — composition**: the boundary laws of a product equal the sums of the operands' laws.
- **Property**: `Empty` arises only from disjoint hard domains.

## Code references

- `packages/treetime/src/timetree/inference/branch_length_likelihood.rs:62-81` — `create_simple_grid`
- `packages/treetime/src/timetree/inference/backward_pass.rs:68-80, 110-122` — child fold, message construction
- `packages/treetime/src/timetree/inference/forward_pass.rs:98-130` — cavity division, forward message
- `packages/treetime/src/timetree/convergence/likelihood.rs:29-70` — `compute_positional_log_lh`
- `packages/treetime-distribution/src/distribution_ops/multiply.rs:138-157` — `multiply_function_function`
- `packages/treetime-distribution/src/distribution_ops/divide.rs` `fn divide_function_by_function()` -- tail-aware bounds (now shares the multiplication support rule)
- `packages/treetime-distribution/src/distribution_ops/time_bounds.rs:52-64` — `distribution_support_intersection`
- `packages/treetime-distribution/src/policy.rs` — `YAxisPolicy`, `supports_zero_boundary`, `SupportsConvolution`
- `packages/treetime-distribution/src/distribution_core/distribution.rs:129-133, 214-243` — `time_bounds` (panics on `Empty`), `max_value` / `scale_by` / `normalize`
- `packages/treetime-distribution/src/distribution_core/function.rs:297-306` — `likely_time` via `argmax`
- `packages/treetime-grid/src/grid_fn.rs:25-36, 67-80, 485-497` — `BoundaryBehavior`, policy reset in `from_grid_array`, policy preservation in `resample`
- `packages/legacy/treetime/treetime/distribution.py:216` — v0 neg-log storage
- `packages/legacy/treetime/treetime/node_interpolator.py:231-256` — v0 slope-based convolution tails
- `packages/legacy/treetime/treetime/branch_len_interpolator.py:36-60` — v0 grid reaching `MAX_BRANCH_LENGTH`
