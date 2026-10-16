# Brent per-edge optimizer ignores its termination status and uses an absolute tolerance

> [!IMPORTANT]
> **Decision required.** Two choices remain: what to do when Brent stops at its iteration limit (today the result is used silently), and whether the stopping tolerance stays absolute in $s = \sqrt{t}$ (argmin default) or becomes relative as in v0. The measured effects on the six example datasets are about $10^{-6}$ log-likelihood units per edge, so the choice is about guarantees, not about current outputs.

## Problem

`fn brent_sqrt_inner()` [packages/treetime/src/optimize/method_brent.rs#L12-L35](../../packages/treetime/src/optimize/method_brent.rs#L12-L35) runs `argmin::solver::brent::BrentOpt` on $-\ell(s^2)$ over the bracket of `fn brent_bracket()` [method_brent.rs#L37-L42](../../packages/treetime/src/optimize/method_brent.rs#L37-L42) with `max_iters(BRENT_MAX_ITER = 50)` [method_brent.rs#L10](../../packages/treetime/src/optimize/method_brent.rs#L10), and returns `best_param` [method_brent.rs#L29](../../packages/treetime/src/optimize/method_brent.rs#L29).

- **Termination status is not read.** When the solver stops at `max_iters`, argmin terminates with `TerminationReason::MaxItersReached`, and the result is used as if it had converged. No warning, no fallback
- **The stopping tolerance is absolute in $s$.** `BrentOpt::new()` sets `eps = sqrt(f64::EPSILON)` ($1.49\cdot10^{-8}$) and `t = 1e-5`, and stops when the bracket is within $\text{tol} = \text{eps} \cdot |s| + 10^{-5}$ of the current point (argmin revision `f5c3fd3`, `crates/argmin/src/solver/brent/brentopt.rs`, `fn new()` and `fn next_iter()`). For short branches the absolute term dominates: at $t = 10^{-4}$ ($s = 10^{-2}$) the accuracy is about $3\cdot10^{-5}$ in $s$, or 0.6% of $t$
- **The result can be worse than the start.** The start only sets the bracket upper end `1.5 t + one_mutation`; argmin never evaluates it. The first sample is the golden-section point $a + c (b - a)$ of the bracket, and the solver returns the best point it sampled. When the start is already within the tolerance of the maximum, the result can have a slightly lower $\ell$

## Evidence

Throwaway instrumentation of `fn optimize_edge()` on `rust` at `77dd3a84`, with the reference maximum described in [M-optimize-newton-returns-unverified-point.md](M-optimize-newton-returns-unverified-point.md#evidence).

Measured: `treetime optimize` with default settings on `flu/h3n2/{20,200}`, `ebola/{20,100}`, `zika/86`, `rsv/a/100`, for `--dense true|false` and `--model infer|jc69` (24 Brent runs, 9,542 per-edge optimizations):

- **Termination**: every solver run ended with `SolverConverged`; the largest iteration count was 29 (mean 12.3), so the limit of 50 was never reached
- **Final $\ell$ below the start**: 875 edges (9%), by at most $9.6\cdot10^{-7}$. On no edge does the start exceed the reference maximum by more than $10^{-9}$, so the reference bounds the start and the two figures below are consistent
- **Final $\ell$ below the reference maximum**: 61 edges by more than $10^{-6}$, largest shortfall $1.1\cdot10^{-6}$
- **Near the bracket lower end**: 752 edges returned a point within $10^{-4}$ in $s$ of the lower end $s = 10^{-6}$ (the probe's criterion); on all of them `reconcile_zero_boundary` replaced it with zero (see [N-optimize-zero-boundary-handled-in-several-places.md](N-optimize-zero-boundary-handled-in-several-places.md)). At full precision these points lie at $t$ between $1.5\cdot10^{-10}$ and $4.4\cdot10^{-10}$, so $s$ between $1.2\cdot10^{-5}$ and $2.1\cdot10^{-5}$: inside the solver's absolute tolerance $10^{-5}$ of the endpoint, not on it, because the solver never samples its endpoints
- **Multimodality**: no edge had more than one local maximum on the reference grid

> [!IMPORTANT]
> **Investigation required.** Not measured: `mpox/clade-ii/100` (run stopped partway), `dengue/100`, datasets above 200 tips, long genomes where $\text{one\_mutation}$ is below $10^{-5}$ (for example `mpox`, `tb`, `sc2`) and the absolute tolerance is coarser than one mutation, fixed non-unimodal models (`k80`, `hky85`, `tn93`), and the `timetree` command. The partial `mpox/clade-ii/100` run showed a shortfall of $1.1\cdot10^{-5}$ on one edge with $t \approx 3.9\cdot10^{-6}$, which suggests a larger shortfall on long genomes, where the absolute tolerance is a larger fraction of $t$.

## Comparison with v0

- `GTR.optimal_t_compressed()` [packages/legacy/treetime/treetime/gtr.py#L816-L920](../../packages/legacy/treetime/treetime/gtr.py#L816-L920) calls `scipy.optimize.minimize_scalar(method='brent')` with the three-point bracket $(-\sqrt{4}, \sqrt{h}, \sqrt{4})$ in $s$, where $h$ is the profile Hamming distance and $4$ is `MAX_BRANCH_LENGTH`. SciPy takes the outer points as limits of the search (every step stays inside $[-2, 2]$, so $t \leq 4$), starts at the middle point, and raises an error unless the objective at the middle point is below both outer points. The result is never worse than the start $\sqrt{h}$. On the profile path the objective adds the penalty $e^{s^4 / 10^4}$ [gtr.py#L867](../../packages/legacy/treetime/treetime/gtr.py#L867); the compressed path has none. $s = 0$ is an interior point
- The tolerance is relative: `optimize_tree_marginal()` passes `tol = 1e-8 + 0.01 ** (i + 1)` for iteration `i` [packages/legacy/treetime/treetime/treeanc.py#L1313](../../packages/legacy/treetime/treetime/treeanc.py#L1313) ($0.01000001$, then $0.00010001$, then $0.00000101$, ...), and scipy Brent uses $\text{tol} \cdot |s| + 10^{-11}$
- v0 reads the termination status: when `opt['success']` is false it returns the Hamming distance instead [gtr.py#L916-L918](../../packages/legacy/treetime/treetime/gtr.py#L916-L918). scipy's iteration limit is 500 (`maxiter=500`)

## Options

1. **Read the termination status**: on `MaxItersReached`, raise a `RunWarning`; for the substitution-only path, v0's fallback is the Hamming distance. Keeping the better of the start and the result differs from v0 and needs approval. Does not change outputs on the measured datasets
2. **Relative tolerance**: call `BrentOpt::set_tolerance(eps, t)` with a small absolute part (for example $10^{-11}$, as scipy) so the accuracy scales with $s$. Changes outputs at the $10^{-6}$ level; needs a decision whether to also copy the v0 per-iteration schedule (see [N-optimize-pre-step-coarse-tolerance.md](N-optimize-pre-step-coarse-tolerance.md)). A tolerance change alone does not restore v0 parity, because the objective (penalty) and the bounds ($t \leq 4$ against `1.5 t + one_mutation`) also differ
3. **Accept-if-not-worse**: return the start when $\ell(\text{result}) < \ell(\text{start})$. Removes the 875 small regressions; costs one evaluation per edge. v0 reaches the same property through its start point, but v0 starts from the Hamming distance, not from the current branch length

## Cross-references

- [N-optimize-per-edge-constants-without-provenance.md](N-optimize-per-edge-constants-without-provenance.md): `BRENT_MAX_ITER`, bracket constants
- [N-optimize-pre-step-coarse-tolerance.md](N-optimize-pre-step-coarse-tolerance.md)
