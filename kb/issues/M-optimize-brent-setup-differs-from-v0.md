# Brent per-edge optimizer setup differs from v0

> [!IMPORTANT]
> **Research and decision required before any change.** The v1 Brent setup stays as it is until the maintainers decide. This issue records how the v1 and v0 setups differ, axis by axis, and which setup is better on each axis with the scientific reason. The open choices are in [Options](#options). The measured effect on the six example datasets is below $10^{-6}$ log-likelihood units per edge; larger effects are possible on long genomes, long branches, and models that were not measured (see the **Investigation required** blocks).

Both versions optimize each branch length with Brent's method in $s = \sqrt{t}$, but they differ in the search interval, the start, the objective, the stopping rule, and the failure handling. `--opt-method brent` (the default) therefore does not reproduce v0 per-edge results exactly.

Notation: $t$ is the branch length, $s = \sqrt{t}$, $\ell$ the log-likelihood of the edge, $L$ the number of alignment sites, $h$ the Hamming distance between the parent and child profiles, $t_0$ the branch length before the optimization step.

## Locations

- v1 solver: `fn solve_brent()` [packages/treetime/src/optimize/method_brent.rs#L17-L28](../../packages/treetime/src/optimize/method_brent.rs#L17-L28), termination status `fn brent_sqrt_minimize()` [method_brent.rs#L30-L63](../../packages/treetime/src/optimize/method_brent.rs#L30-L63), upper-end check `fn ascends_at_upper_end()` [method_brent.rs#L65-L79](../../packages/treetime/src/optimize/method_brent.rs#L65-L79), bounds `fn brent_bracket()` [method_brent.rs#L81-L86](../../packages/treetime/src/optimize/method_brent.rs#L81-L86), `BRENT_MAX_ITER` [method_brent.rs#L11](../../packages/treetime/src/optimize/method_brent.rs#L11), `BRENT_REL_TOL` and `BRENT_ABS_TOL` [method_brent.rs#L13-L15](../../packages/treetime/src/optimize/method_brent.rs#L13-L15)
- v1 caller: `fn optimize_edge()` [packages/treetime/src/optimize/dispatch.rs#L225-L261](../../packages/treetime/src/optimize/dispatch.rs#L225-L261), start `fn start_branch_length()` [packages/treetime/src/optimize/zero_boundary.rs#L13-L30](../../packages/treetime/src/optimize/zero_boundary.rs#L13-L30), lower bound `fn min_branch_length_for_indels()` [zero_boundary.rs#L32-L38](../../packages/treetime/src/optimize/zero_boundary.rs#L32-L38), upper bound `fn search_upper_bound()` [zero_boundary.rs#L120-L122](../../packages/treetime/src/optimize/zero_boundary.rs#L120-L122) with `MIN_SEARCH_UPPER_BOUND = 0.5` [zero_boundary.rs#L11](../../packages/treetime/src/optimize/zero_boundary.rs#L11), non-convergence warning `ConvergenceSummary::report()` [packages/treetime/src/optimize/solver_outcome.rs#L43](../../packages/treetime/src/optimize/solver_outcome.rs#L43)
- argmin `BrentOpt` 0.11.0 (the workspace pin `=0.11.0` from crates.io, tag `argmin-v0.11.0`): constructor defaults [brentopt.rs#L66-L71](https://github.com/argmin-rs/argmin/blob/c94c32adefd6c2525ce05806092ca868ec85fba4/crates/argmin/src/solver/brent/brentopt.rs#L66-L71), first sample [brentopt.rs#L109-L123](https://github.com/argmin-rs/argmin/blob/c94c32adefd6c2525ce05806092ca868ec85fba4/crates/argmin/src/solver/brent/brentopt.rs#L109-L123), stopping test [brentopt.rs#L133-L143](https://github.com/argmin-rs/argmin/blob/c94c32adefd6c2525ce05806092ca868ec85fba4/crates/argmin/src/solver/brent/brentopt.rs#L133-L143)
- v0 solver: `GTR.optimal_t_compressed()` [packages/legacy/treetime/treetime/gtr.py#L816-L920](../../packages/legacy/treetime/treetime/gtr.py#L816-L920), penalty [gtr.py#L867](../../packages/legacy/treetime/treetime/gtr.py#L867), fallback [gtr.py#L916-L918](../../packages/legacy/treetime/treetime/gtr.py#L916-L918), `MAX_BRANCH_LENGTH = 4.0` [packages/legacy/treetime/treetime/config.py#L34-L36](../../packages/legacy/treetime/treetime/config.py#L34-L36)
- v0 caller: `TreeAnc.optimize_tree_marginal()` tolerance schedule [packages/legacy/treetime/treetime/treeanc.py#L1313](../../packages/legacy/treetime/treetime/treeanc.py#L1313), profile calls [treeanc.py#L1327-L1341](../../packages/legacy/treetime/treetime/treeanc.py#L1327-L1341)
- SciPy 1.18.1 (the version in `pixi.toml`): `minimize_scalar` passes `tol` as `xtol` to Brent [_minimize.py#L1011](https://github.com/scipy/scipy/blob/v1.18.1/scipy/optimize/_minimize.py#L1011); three-point bracket check [_optimize.py#L2471-L2489](https://github.com/scipy/scipy/blob/v1.18.1/scipy/optimize/_optimize.py#L2471-L2489); start and stopping test [_optimize.py#L2503-L2526](https://github.com/scipy/scipy/blob/v1.18.1/scipy/optimize/_optimize.py#L2503-L2526); `_mintol = 1e-11` [_optimize.py#L2447](https://github.com/scipy/scipy/blob/v1.18.1/scipy/optimize/_optimize.py#L2447); success flag [_optimize.py#L2738](https://github.com/scipy/scipy/blob/v1.18.1/scipy/optimize/_optimize.py#L2738)

## Summary

| Axis                          | v0                                                                                          | v1                                                                                                       | Better                                                           |
| ----------------------------- | ------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------- |
| Search interval and start     | $s \in [-2, 2]$ ($t \le 4$), start $\sqrt{h}$, result never worse than start                | $t \in [\max(t_\text{min}, 10^{-12}),\ \max(1.5\,t_0 + \text{one mutation}, 0.5)]$, golden-section start | Start: v0. Upper bound: neither. Zero: equal outcome, v0 simpler |
| Objective                     | $-\ell(s^2) + e^{s^4/10^4}$ on the profile path                                             | $-\ell(s^2)$                                                                                             | v1                                                               |
| Tolerance and iteration limit | $\text{tol}_i \cdot \lvert s \rvert + 10^{-11}$, $\text{tol}_i = 10^{-8} + 0.01^{i+1}$; 500 | $\sqrt{\epsilon}\,\lvert s \rvert + 10^{-5}$; 50, status read                                            | Tolerance shape: v1. Status: equal                               |
| Failure handling              | Non-convergence: Hamming distance. Invalid bracket: uncaught `ValueError`                   | Cost-function error: `Err` to the caller. Iteration limit or upper-end stop: best point kept, run warning | Errors: v1. Non-convergence: v1                                  |

## Axis 1: search interval and start

### What each version does

- **v0** passes the three-point bracket $(-\sqrt{4},\ \sqrt{h},\ \sqrt{4})$ in $s$. SciPy treats the outer points as the limits of the search: every step stays inside $[-2, 2]$, so $t \le 4$. The objective is even in $s$, so the negative half repeats the positive half, and $s = 0$ ($t = 0$) is an interior point. Brent starts at the middle point $\sqrt{h}$, where $h$ is the profile Hamming distance ($h \le 1$), and SciPy raises a `ValueError` unless $f(\sqrt{h}) < f(\pm 2)$. Brent keeps the best point it has evaluated, and the start is evaluated, so the result is never worse than $\sqrt{h}$
- **v1** passes $[\sqrt{\text{lower}}, \sqrt{\text{upper}}]$ with $\text{lower} = \max(t_\text{min}, 10^{-12})$, where $t_\text{min} = 0.01 \cdot \text{one mutation}$ on edges whose indel term needs a positive length and $0$ otherwise, and $\text{upper} = \max(1.5\,t_0 + \text{one mutation}, 0.5)$. argmin takes no start: its first sample is the golden-section point $a + 0.382\,(b - a)$ of the interval, and $t_0$ only sets the upper end. Zero is outside the interval; `fn optimize_edge()` handles it before the solver (`fn is_zero_branch_optimal()`) and after it (`fn reconcile_zero_boundary()`), see [N-optimize-zero-boundary-handled-in-several-places.md](N-optimize-zero-boundary-handled-in-several-places.md)

### Consequences

- **Result below the start (v1 only).** When $t_0$ is already within the tolerance of the maximum, v1 can return a point with a lower $\ell$ than $t_0$. Measured on 875 of 9,542 edges (9%), by at most $9.6\cdot10^{-7}$ (see [Evidence](#evidence))
- **Long branches.** v0 cannot return $t > 4$; a saturated branch stops at 4 (v0 logs a warning above $0.9 \cdot 4$). v1 has no fixed cap, but one step can reach at most $1.5\,t_0 + \text{one mutation}$ (at least 0.5), so a branch whose maximum lies far above $t_0$ needs several outer iterations to reach it: from $t_0 = 0.01$ to a maximum at $t = 2$ needs the upper ends $0.5, 0.75, 1.1, 1.7, 2.5$ (ignoring the one-mutation term and damping). v1 reports such a stop: when the converged point lies within $3\,\text{tol}$ of $\sqrt{\text{upper}}$ and the log-likelihood still increases at the end (`fn ascends_at_upper_end()`), the solution gets the status `Convergence::SearchBound`, and the `branch-lengths-not-converged` run warning names the edge
- **Start outside the interval.** Brent evaluates only points inside $[\sqrt{\text{lower}}, \sqrt{\text{upper}}]$, so it never evaluates a start below its lower bound. The safeguarded Newton method evaluates a start outside the interval, keeps it as the incumbent when its log-likelihood is finite, and returns it when the search finds no better point (`fn feasible_start_outside()` [packages/treetime/src/optimize/method_newton.rs#L64-L75](../../packages/treetime/src/optimize/method_newton.rs#L64-L75))
- **Near-zero branches.** Both reach $t = 0$ on the measured data: v0 because 0 is interior, v1 through the zero handling around the solver (752 measured edges ended within $10^{-4}$ in $s$ of the lower end and were replaced by zero)

### Established tools

- **IQ-TREE** (revision `63c330d`): Newton-Raphson is the default; with `-brent`, `PhyloTree::optimizeOneBranch()` calls `Optimization::minimizeOneDimen()` in $t$ (not $\sqrt{t}$) on $[t_\text{min}, 10]$, where $t_\text{min} = 10^{-6}$, or $0.1 / L$ for alignments of at least $10^5$ sites when the user sets no minimum. It starts at the current length and returns the start when the Brent result is worse ([optimization.cpp#L322-L373](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/optimization.cpp#L322-L373), [phylotree.cpp#L2669-L2675](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/tree/phylotree.cpp#L2669-L2675), [phyloanalysis.cpp#L5259-L5262](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/main/phyloanalysis.cpp#L5259-L5262), [tools.cpp#L7304-L7307](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/tools.cpp#L7304-L7307))
- **RAxML-NG** (revision `d396351`) with coraxlib optimizes branch lengths with Newton-Raphson (default `CORAX_OPT_BLO_NEWTON_FAST`) on $[10^{-6}, 100]$; it uses Brent only for one-dimensional model parameters ([constants.hpp#L29-L31](https://github.com/amkozlov/raxml-ng/blob/d396351ee1263b704adb135edd6a4a84522dbd5b/src/constants.hpp#L29-L31), [Options.cpp#L39](https://github.com/amkozlov/raxml-ng/blob/d396351ee1263b704adb135edd6a4a84522dbd5b/src/Options.cpp#L39))
- PAML and PhyML were not consulted

### Which is better

- **Start: v0.** A line search that starts at a known point and never returns a worse point gives a monotone non-decreasing likelihood per edge update. IQ-TREE has the same property through its explicit rollback. v1 lacks it. v1's choice of $t_0$ as the reference (rather than v0's $h$) is the better start, because $t_0$ is usually closer to the maximum after the first outer iteration
- **Upper bound: neither.** v0's cap of 4 is below the caps of IQ-TREE (10) and RAxML-NG (100) and truncates saturated branches. v1's adaptive bound has no cap but delays long branches by several outer iterations
- **Zero: equal outcome on the measured data, v0 simpler.** v0 reaches zero inside one solver call; v1 needs three mechanisms around the solver

> [!IMPORTANT]
> **Investigation required.** The `SearchBound` status makes upper-bound stops visible in the `branch-lengths-not-converged` run warning, but their frequency on the example datasets is not measured. Also not measured: whether the outer loop stops before such a branch reaches its maximum, and how many v0 branches stop at the cap $t = 4$ on the example datasets.

## Axis 2: objective

### What each version does

- **v0** minimizes $-\ell(s^2) + e^{s^4/10^4}$ when it optimizes with profiles ([gtr.py#L867](../../packages/legacy/treetime/treetime/gtr.py#L867)); `optimize_tree_marginal()` always uses the profile path. In the code the argument `t` of `_neg_prob` is $s$, so the penalty is $e^{t^2/10^4}$ in branch-length units. The path with fixed sequences (`profiles=False`) has no penalty
- **v1** minimizes $-\ell(s^2)$, including the Poisson indel term, without a penalty

### Consequences

- **Size of the shift.** With $e^{x} \approx 1 + x$, the penalty adds the gradient $4s^3/10^4$ in $s$. With the curvature $4L$ of $\ell$ in $s$ (derivation in Axis 3), the maximum moves by about $\delta s \approx s^3 / (10^4 L)$: $10^{-7}$ at $t = 1$ and $L = 10^3$. The shift matters only for long branches on short alignments
- **Bracket validity.** The penalty makes $f(\pm 2)$ larger than $f(\sqrt{h})$ when $\ell$ is flat (an edge without informative sites), so SciPy's bracket condition holds there. v1 needs no such condition, because argmin's `BrentOpt` takes an interval, not a three-point bracket

### Which is better

**v1.** The quantity estimated is the maximum-likelihood branch length; a penalty makes it a penalized estimate whose value v0 then reports without the penalty. IQ-TREE and RAxML-NG optimize the unpenalized likelihood inside fixed bounds. v1 needs a finite upper bound to keep saturated branches finite, and it has one (Axis 1).

> [!IMPORTANT]
> **Investigation required.** The effect of the penalty on v0 outputs was not measured; the estimate above is a derivation under the Poisson approximation.

## Axis 3: tolerance, iteration limit, and termination status

### What each version does

- **v0**: SciPy's Brent stops when the bracket is within $\text{tol}_1 = \text{tol} \cdot |s| + 10^{-11}$ of the current point, with $\text{tol}_i = 10^{-8} + 0.01^{i+1}$ in outer iteration $i$: $0.01000001$, then $0.00010001$, then $0.00000101$, and so on. The timetree pre-step of v0 runs one outer iteration (`optimize_tree(max_iter=1)` [packages/legacy/treetime/treetime/treetime.py#L243](../../packages/legacy/treetime/treetime/treetime.py#L243)), so it uses only $\text{tol}_0 \approx 0.01$. The limit is 500 iterations; `success` is false when the limit is reached or the result is NaN, and v0 reads it (Axis 4)
- **v1**: `BRENT_REL_TOL` $= \sqrt{\epsilon} \approx 1.49\cdot10^{-8}$ and `BRENT_ABS_TOL` $= 10^{-5}$, passed to `BrentOpt::set_tolerance()` and equal to argmin's defaults, give $\text{tol} = \sqrt{\epsilon}\,|s| + 10^{-5}$, the same in every outer iteration and in the timetree pre-step, because `fn optimize_edge_lengths()` takes no tolerance [packages/treetime/src/optimize/dispatch.rs#L22-L29](../../packages/treetime/src/optimize/dispatch.rs#L22-L29). Since $|s| \le \sqrt{\text{upper}}$, the relative part is at most $1.5\cdot10^{-7}$ for every branch up to $t = 100$, two orders of magnitude below $10^{-5}$, so the tolerance is in practice absolute: $10^{-5}$ in $s$. The limit is `BRENT_MAX_ITER = 50`; `fn brent_sqrt_minimize()` reads the termination reason, keeps the best point on `MaxItersReached`, and the optimization run raises one `branch-lengths-not-converged` run warning that names the affected edges (together with the edges that stopped at the upper end, Axis 1). The returned point is within about $3\,\text{tol}$ of the local minimum (argmin documentation of `set_tolerance()`)

### Accuracy relative to the statistical error

For a short branch with $k$ differences over $L$ sites, the Poisson approximation gives $\ell(t) \approx k \ln t - L t$, so in $s$:

$$
\ell(s) \approx 2k \ln s - L s^2, \qquad \frac{d^2 \ell}{ds^2} = -\frac{2k}{s^2} - 2L = -4L \ \text{ at } \hat{s}^2 = k/L
$$

The standard error of $\hat{s}$ is therefore $1/(2\sqrt{L})$, independent of $t$: the square-root transform stabilizes the variance of a Poisson count. Two consequences follow:

- **v1 (absolute in $s$).** An error $\delta s \le 3\cdot10^{-5}$ is $6\cdot10^{-5}\sqrt{L}$ standard errors for every $t$: 0.6% at $L = 10^4$, 6% at $L = 10^6$, 13% at $L = 4.4\cdot10^6$. The log-likelihood shortfall is at most $\tfrac{1}{2} \cdot 4L \cdot \delta s^2 \approx 1.8\cdot10^{-9} L$. The measured largest shortfall below the reference maximum, $1.1\cdot10^{-6}$ on alignments of at most about $2\cdot10^4$ sites, and the $1.1\cdot10^{-5}$ of the partial `mpox/clade-ii/100` run (about $2\cdot10^5$ sites) are below this bound
- **v0 (relative in $s$).** An error $\delta s \le 3\,\text{tol}_i\,s$ is $6\,\text{tol}_i\,s\sqrt{L}$ standard errors. In the first outer iteration ($\text{tol}_0 \approx 0.01$) this is $0.6$ standard errors at $t = 0.01$, $L = 10^4$, which damping then mixes with the old length; from the third iteration on, v0 is more accurate than v1 at every $t \le 4$. Near $t = 0$ the relative tolerance resolves $s$ to $10^{-11}$, far below any statistical meaning

The Poisson approximation ignores multiple hits and ambiguous sites, which lower the curvature; the figures are orders of magnitude, not bounds for every model.

### Established tools

- **IQ-TREE** Brent: $\text{tol}_1 = t_\text{min} \cdot |t| + 10^{-10}$ in $t$, at most 100 iterations; on the limit it returns the best point without a report ([optimization.cpp#L174-L215](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/optimization.cpp#L174-L215))
- **RAxML-NG**: Newton-Raphson with tolerance $10^{-7}$ in $t$ (`RAXML_BRLEN_TOLERANCE`)

### Which is better

- **Tolerance shape: v1.** An absolute tolerance in $s$ is a constant fraction of the standard error at every branch length, which is the right scale for a statistical estimate. v0's relative tolerance is coarse in the first outer iteration and over-resolves near zero. Neither scales the tolerance with $L$: an absolute tolerance proportional to $1/\sqrt{L}$ in $s$ would hold the error at a fixed fraction of the standard error for every alignment length
- **Iteration limit: equal in practice.** Neither limit was reached on the measured data (v1 largest count 29 of 50)
- **Termination status: equal.** Both read it; they differ in what they do at the limit (Axis 4)

## Axis 4: failure handling

### What each version does

- **v0**: when `opt['success']` is false (500 iterations or NaN), it returns the Hamming distance $h$ ([gtr.py#L916-L918](../../packages/legacy/treetime/treetime/gtr.py#L916-L918)). When the bracket condition fails, SciPy raises `ValueError`; `optimal_t_compressed()` catches only `ImportError`, so the exception ends the run
- **v1**: an error of the cost function or of argmin becomes an `Err` with the branch length in the message and ends the run; a stop at `MaxItersReached` or at the upper end with a rising log-likelihood (`SearchBound`, Axis 1) returns the best point and raises a run warning (Axis 3)

### Consequences

- $h$ is a model-free proportion of differing sites. It ignores multiple substitutions and underestimates $t$ (under JC69, $t = -\tfrac{3}{4}\ln(1 - \tfrac{4}{3}h) > h$), and it discards the best point Brent found. The substitution is silent apart from a log line at verbosity 4

### Which is better

- **Errors: v1.** A visible error is better than v0's silent Hamming substitution, and v0's uncaught `ValueError` is an error too, without context
- **Non-convergence: v1.** The best point Brent found is the better estimate (IQ-TREE keeps it), and the run warning tells the user which edges did not converge

## Evidence

Measurements were taken with throwaway instrumentation of `fn optimize_edge()` on `rust` at `77dd3a84`. For every edge and optimize iteration, the probe recorded the exit, the start, the solver output, the final output after `fn reconcile_zero_boundary()`, and a reference maximum: $\ell$ on 800 log-spaced points over $[10^{-4} \cdot \text{one\_mutation}, 10]$ plus $t = 0$, with golden-section refinement around each grid local maximum. A finite reference grid does not prove a global maximum.

Measured: `treetime optimize` with default settings on `flu/h3n2/{20,200}`, `ebola/{20,100}`, `zika/86`, `rsv/a/100`, for `--representation dense|sparse` and `--model infer|jc69` (24 Brent runs, 9,542 per-edge optimizations):

- **Termination**: every solver run ended with `SolverConverged`; the largest iteration count was 29 (mean 12.3), so the limit of 50 was never reached. The probe ran before the upper-end check existed, so it does not separate converged points from stops at the upper end
- **Final $\ell$ below the start**: 875 edges (9%), by at most $9.6\cdot10^{-7}$. On no edge does the start exceed the reference maximum by more than $10^{-9}$, so the reference bounds the start and the two figures below are consistent
- **Final $\ell$ below the reference maximum**: 61 edges by more than $10^{-6}$, largest shortfall $1.1\cdot10^{-6}$
- **Near the interval lower end**: 752 edges returned a point within $10^{-4}$ in $s$ of the lower end $s = 10^{-6}$ (the probe's criterion); on all of them `reconcile_zero_boundary` replaced it with zero. At full precision these points lie at $t$ between $1.5\cdot10^{-10}$ and $4.4\cdot10^{-10}$, so $s$ between $1.2\cdot10^{-5}$ and $2.1\cdot10^{-5}$: inside the solver's absolute tolerance $10^{-5}$ of the endpoint, not on it, because the solver never samples its endpoints
- **Multimodality**: no edge had more than one local maximum on the reference grid

> [!IMPORTANT]
> **Investigation required.** Not measured: `mpox/clade-ii/100` (run stopped partway), `dengue/100`, datasets above 200 tips, long genomes where the absolute tolerance is a larger fraction of the standard error (for example `mpox`, `tb`, `sc2`), fixed non-unimodal models (`k80`, `hky85`, `tn93`), and the `timetree` command. The partial `mpox/clade-ii/100` run showed a shortfall of $1.1\cdot10^{-5}$ on one edge with $t \approx 3.9\cdot10^{-6}$. No v0 run was instrumented, so the per-axis effect on the difference between v0 and v1 outputs is not measured.

## Possible downstream effect

[M-optimize-gm-per-branch-divergence.md](M-optimize-gm-per-branch-divergence.md) records per-branch differences from the v0 golden master beyond 10%. The setup differences above may contribute, but no measurement links them; the measured per-edge effects here are far smaller than that divergence.

## Options

Each option is independent of the others. v1 reads the termination status (Axis 3) and keeps the best point at the iteration limit instead of v0's fallback to $h$, which is scientifically worse (Axis 4).

1. **Accept-if-not-worse**: evaluate $\ell(t_0)$ and return $t_0$ when $\ell(\text{result}) < \ell(t_0)$, as IQ-TREE does. Removes the 875 small regressions; costs one evaluation per edge. v0 reaches the same property through its start point, but v0 starts from $h$, not from $t_0$
2. **Tolerance**: keep the absolute tolerance in $s$ and scale it with the alignment, for example $c / (2\sqrt{L})$ in place of the constant `BRENT_ABS_TOL`; or copy v0's relative tolerance ($10^{-11}$ absolute part, as SciPy) and its per-iteration schedule, including the coarse first iteration that the timetree pre-step uses. Either changes outputs at about the $10^{-6}$ level on the measured datasets. A tolerance change alone does not restore v0 parity, because the objective and the interval also differ
3. **Upper bound**: a fixed cap (v0: 4, IQ-TREE: 10, RAxML-NG: 100) instead of, or in addition to, the adaptive $1.5\,t_0 + \text{one mutation}$; depends on the investigation in Axis 1
4. **Penalty**: keep v1 without the v0 penalty (recommended by Axis 2), or add it for parity

> [!IMPORTANT]
> **Decision required.** Which of options 1 to 4 to adopt. Option 1 does not move v1 away from v0 on the measured data and closes the guarantee that v0 and IQ-TREE give. Options 2 to 4 change outputs and each is a choice between v0 parity and the better-justified setup of its axis.

## Cross-references

- [N-optimize-per-edge-constants-without-provenance.md](N-optimize-per-edge-constants-without-provenance.md): `BRENT_MAX_ITER`, interval constants
- [N-optimize-zero-boundary-handled-in-several-places.md](N-optimize-zero-boundary-handled-in-several-places.md)
- [N-optimize-grid-search-upper-bound-capped.md](N-optimize-grid-search-upper-bound-capped.md)
- [kb/decisions/optimize-newton-raphson-per-edge.md](../decisions/optimize-newton-raphson-per-edge.md)
