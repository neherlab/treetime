# Branch Length Optimization

v0 has no standalone optimize command. `TreeAnc.optimize_tree()` and `optimize_tree_marginal()` perform branch length optimization inline during timetree and ancestral workflows. v1 extracts this into a standalone `optimize` command and a reusable `PartitionOptimizeOps` trait system shared with timetree.

## Standalone Command (v1-only)

- [x] `optimize` command with `--tree`, `--aln`, `--outdir`
- [x] Mixed dense and sparse partition setup (one of each, from same alignment)
- [x] Initial branch-length guess from observed mutation counts (excludes deletion and ambiguous positions)
- [x] `--branch-length-initial-guess` flag: `auto` (selective fill, treats zero BL as invalid when indels present), `always` (overwrite all), `never` (error on missing/NaN, and on zero branch length when indels are present)
- [x] Output annotated Newick and Nexus trees
- [x] GTR parameters written to JSON

## Per-Edge Likelihood (v0 parity via different method)

v0 uses Brent's method (`scipy.optimize.minimize_scalar`) in sqrt(t) space with Hamming distance bracket. v1 offers six methods via `--opt-method`: Newton and Brent in three parameterizations ($t$, $\sqrt{t}$, $\ln(t)$). Default is `brent-sqrt` (matches v0). Both v0 and v1 use eigendecomposition-based likelihood (`expQt = V diag(exp(lambda*t)) V_inv`).

- [x] Eigenvalue-space coefficient caching (dense: `msg.dot(V) * msg.dot(V_inv.T)`, sparse: per-site with multiplicity)
- [x] Analytical first and second derivatives (v1-only, v0 uses derivative-free Brent)
- [x] Newton's method with clamped step (max 10 inner iterations, step in `[-1.0, bl]`, absolute tolerance floor 1e-8 subs/site)
- [x] Newton's method in sqrt(t) space (`newton_sqrt_inner`, reduces indel Hessian singularity from O(1/t^2) to O(1/t))
- [x] Newton's method in ln(t) space (`newton_log_inner`, eliminates indel Hessian singularity entirely, bounded curvature -mu\*t)
- [x] Brent's method in t space (`brent_inner`, derivative-free via `argmin::BrentOpt`)
- [x] Brent's method in sqrt(t) space (`brent_sqrt_inner`, default, matches v0 exactly)
- [x] Brent's method in ln(t) space (`brent_log_inner`, smoothest objective surface)
- [x] Per-edge optimization method selection (`--opt-method`: 6 methods, default `brent-sqrt`)
- [x] Grid search fallback when second derivative >= 0 (100 points, log-spaced grid with 0.5 subs/site minimum upper bound)
- [x] Zero branch length short-circuit (combined likelihood > 0.01 and derivative < 0 at zero)
- [x] `compute_derivatives` flag to skip derivative computation for log-likelihood-only evaluation
- [x] Collect dense contribution (`PartitionMarginalDense`)
- [x] Collect sparse contribution (`PartitionMarginalSparse`, multiplicity-weighted)
- [x] Unified mixed-partition evaluation (`evaluate_mixed()` sums metrics across partition types)
- [x] sqrt(t) reparameterization (Brent in s=sqrt(t) space, matching v0; also ln(t) space variant)
- [ ] Regularization penalty for profile-based optimization (v0: `exp(t^4/10000)` prevents unbounded growth)
- [ ] Hamming distance fallback when optimization fails (v0 falls back to observed distance)

## Convergence Loop

- [x] Iterative marginal reconstruction + optimization loop bounded by `--max-iter`
- [x] Early stop when absolute likelihood change is below `--dp`
- [x] Oscillation detection: stop when 2-step likelihood difference is below `--dp`
- [x] Worsened detection: stop and revert to best branch lengths when likelihood decreases from peak
- [x] Damping in marginal loop (`--damping`, default 0.75, matching v0: `new*(1-d^i) + old*d^i`)
- [x] Damping floor (`DAMPING_FLOOR = 0.01`) prevents fully undamped late iterations
- [ ] Progressive per-iteration tolerance tightening (v0: `tol = 1e-8 + 0.01^(i+1)`, coarse early, tight late)
- [ ] Bifurcating root special handling (v0 optimizes combined root-children length, preserves ratio)
- [ ] Convergence by sequence change count (v0 joint mode: stops when zero nucleotides change)

## GTR Integration

- [x] `--model` flag wired through `get_gtr_sparse()`/`get_gtr_dense()` for all named models and inference
- [ ] GTR inference integrated into optimization loop (v0: `infer_gtr` parameter re-estimates model per iteration)
- [ ] GTR rate optimization (v0: `optimize_gtr_rate()` optimizes overall substitution rate mu)

## Reuse by Timetree

- [x] `PartitionOptimizeOps` trait implemented by both dense and sparse partitions
- [x] `collect_edge_contributions()` gathers contributions for one edge across partition types
- [x] `compute_branch_length_distribution()` evaluates log-likelihood (substitution + Poisson indel) on grid, converts to time-domain distribution
- [x] `evaluate_with_indels_log_lh_only()` for grid evaluation without derivatives (matches `run_optimize_mixed()`'s edge evaluator)

## Branch Length Modes (v0 has 3 modes, v1 implements marginal only)

- [x] Marginal mode (profile-based optimization with full probability vectors)
- [ ] Joint mode (removed in v1, see [intentional change](../decisions/ancestral-joint-reconstruction-removed.md))
- [ ] Input mode (no optimization, Poisson/Gaussian approximation for short/long branches)
- [ ] Auto-detection of branch length mode (v0: defaults to `input` if max branch > 0.1, else `joint`)

## Short Branch Handling

- [x] Zero-branch-length derivative check before optimization
- [ ] Short branch pruning after optimization (v0 joint mode: inside loop; v0 marginal mode: after loop; prunes when `bl < 0.1 * one_mutation` and P(zero) > 0.1)
- [ ] MIN_BRANCH_LENGTH floor for GTR calculations (v0: `1e-3 * one_mutation`)

## Current Limitations

- [/] Command always builds one sparse and one dense partition from the same full alignment
- [ ] `--dense` (parsed but not wired, `infer_dense()` is a stub returning false)
- [ ] Separate dense-only and sparse-only command modes not exposed
- [ ] Standalone `run_optimize_sparse()` zero-branch inconsistency
