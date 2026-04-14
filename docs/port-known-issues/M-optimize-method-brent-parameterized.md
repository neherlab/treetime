# Per-edge optimization: Brent in sqrt(t) and ln(t) spaces

## Problem

The existing `brent_inner` operates in $t$-space. Brent is derivative-free and conditioning-invariant, so it finds the correct optimum regardless of parameterization. But the shape of the objective surface affects convergence speed: steep $\ln(\mu t)$ curvature near $t = 0$ degrades Brent's parabolic interpolation to golden section search in $t$-space. Reparameterizing to $\sqrt{t}$ or $\ln(t)$ smooths the surface and improves convergence.

v0 uses Brent in $\sqrt{t}$ space (`scipy.optimize.minimize_scalar(method='brent')` at `packages/legacy/treetime/treetime/gtr.py:842-892`). The `brent-sqrt` method matches v0 exactly and is the default.

## Background

### Brent's method

`argmin::BrentOpt` minimizes a scalar function within a bracket $[a, b]$ using golden section search + parabolic interpolation. Needs only function evaluations (no derivatives). Convergence order ~1.325.

### v0 implementation details

v0 brackets in $s$-space: `[-sqrt(MAX_BRANCH_LENGTH), sqrt(hamming_distance), sqrt(MAX_BRANCH_LENGTH)]`. The Hamming distance provides the initial bracket midpoint (starting point for parabolic interpolation). On failure, v0 returns Hamming distance as the branch length fallback. v0 also adds a regularization penalty `exp(t^4/10000)` for marginal profile optimization -- v1 does not replicate this; the grid search upper bound and Brent bracket provide a soft cap instead.

v1's approach differs from v0 in two ways: (1) no Hamming distance midpoint (BrentOpt initializes from the bracket), (2) on failure, v1 returns the current branch length rather than Hamming distance. These are acceptable because BrentOpt's internal initialization is adequate for the smooth phylogenetic likelihood surface, and "current branch length" is a better fallback than a crude distance measure.

### Parameterization

**$\sqrt{t}$ space** ($s = \sqrt{t}$): cost function takes $s$, evaluates $-\ell(s^2)$. Bracket endpoints transform as $\sqrt{\text{lower}}$ and $\sqrt{\text{upper}}$. Result squared back to $t$. Tolerance $\epsilon_s$ in $s$-space maps to $t$-space precision $\approx 2s^* \epsilon_s$ -- tighter near zero (desirable for short branches).

**$\ln(t)$ space** ($u = \ln(t)$): cost function takes $u$, evaluates $-\ell(e^u)$. Bracket transforms as $\ln(\text{lower})$ and $\ln(\text{upper})$. Result exponentiated back. Lower bound requires $t > 0$; use $\max(\text{min\_bl}, 10^{-12})$. Tolerance $\epsilon_u$ is a natural relative tolerance ($dt/t \approx du$). Precedent: coalescent Tc optimizer uses `BrentOpt::new(-20.0, 2.0)` in log-space at `packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs:62`.

### Shared infrastructure

Both functions need `evaluate_with_indels_log_lh_only` from the unified module (same as existing `brent_inner`). This function is currently private -- it needs to become `pub(crate)` to be callable from the Brent module.

Each parameterized Brent function needs its own cost function struct implementing `argmin::CostFunction`. The existing `BranchLengthCostFunction` (t-space) is the template -- the only difference is the transform inside `cost()`.

The upper bound computation `max(1.5 * bl + one_mutation, 0.5)` is shared across all three Brent variants.

## Approach

Add two new functions alongside the existing `brent_inner` (which is kept for the `brent` CLI variant). Each follows the same pattern: transform bracket endpoints, create a cost function struct that transforms the parameter inside `cost()`, run BrentOpt, transform result back to $t$-space.

BrentOpt's built-in bracket convergence is sufficient -- no custom tolerance parameters needed for the Brent variants.

## Verification

`./dev/docker/run ./dev/dev t` -- existing tests pass. New functions are not called until [M-optimize-method-dispatch](M-optimize-method-dispatch.md) wires them.

Manual check after dispatch is wired: run `./dev/docker/run ./dev/dev r treetime -- optimize --opt-method=brent-sqrt --tree=data/flu/h3n2/20/tree.nwk --aln=data/flu/h3n2/20/aln.fasta.xz --outdir=tmp/opt-brent-sqrt` and verify completion.

## Dependencies

- Depends on: M-optimize-method-scaffolding (done) -- file split puts Brent code in its own module
- Depended on by: [M-optimize-method-dispatch](M-optimize-method-dispatch.md)

## Cross-references

- v0 implementation: `packages/legacy/treetime/treetime/gtr.py:842-892` (Brent in sqrt-space)
- Coalescent Tc log-space precedent: `packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs:62`
- Current `brent_inner`: `packages/treetime/src/commands/optimize/method_brent.rs`
