# Branch Length Optimization

v0 has no standalone optimize command. `TreeAnc.optimize_tree()` and `optimize_tree_marginal()` perform branch length optimization inline during timetree and ancestral workflows. v1 extracts this into a standalone `optimize` command and a reusable `PartitionOptimizeOps` trait system shared with timetree.

## Standalone Command (v1-only)

- [x] `optimize` command with `--tree`, `--aln`, `--output-all`
- [x] Dense or sparse partition setup, one partition per `--alignment` file, selected by `--dense` for the whole run
- [x] Several alignment files share the branch lengths, with one relative rate per file (`--partition-rates`, default `estimate`: a rate step and a site-weighted normalization after each iteration; `equal` fixes every rate at 1) and one indel rate per file in one combined Poisson term ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md))
- [x] Short-branch floors use the total length of all files
- [x] Tips without a sequence in a file are filled as unknown; the run stops when a file misses more than a third of the tips
- [x] Initial branch-length guess from observed mutation counts (excludes deletion and ambiguous positions)
- [x] `--branch-length-initial-guess` flag: `auto` (selective fill, treats zero BL as invalid when indels present), `always` (overwrite all), `never` (error on missing/NaN, and on zero branch length when indels are present)
- [x] Output annotated Newick and Nexus trees
- [x] Model document (`--output-gtr`) with one entry per alignment file
- [x] Augur-compatible node data JSON (`optimize.augur-node-data.json`, `--output-augur-node-data`): per-node optimized `branch_length` plus `alignment`/`input_tree` metadata, equivalent to `augur refine` without `--timetree`

## Rerooting (v1-only, date-free)

The optimize command has no sampling dates, so only date-free rooting policies are available. Rerooting is opt-in; the input root is kept by default.

- [x] `--reroot=min-dev`: minimize the variance of root-to-tip distances with the clock root search at a fixed rate of zero, every tip at the same date ([kb/algo/reroot.md](../algo/reroot.md)). This implements v0's documented fixed-zero-slope meaning; actual v0 candidate scoring retains an estimated-rate term ([kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md)).
- [x] `--reroot-tips=A,B,...`: root on the branch leading to the MRCA of the named tips (midpoint), no scoring.
- [x] `--keep-root`: explicit form of the default (keep the input root); mutually exclusive with the reroot options.
- [x] Date-dependent methods (`least-squares`, `oldest`, `clock-filter`) rejected at argument validation
- [x] Two-phase pipeline: damped pre-reroot branch-length pass, reroot, partition reconciliation, then re-optimization in the main loop (mirrors the timetree pipeline)

## Per-Edge Likelihood (v0 parity via different method)

v0 uses Brent's method (`scipy.optimize.minimize_scalar`) in sqrt(t) space with Hamming distance bracket. v1 offers two methods via `--opt-method`: Brent and Newton-Raphson, both in $\sqrt{t}$ space. Default is `brent` (matches v0). Both v0 and v1 use eigendecomposition-based likelihood (`expQt = V diag(exp(lambda*t)) V_inv`).

- [x] Eigenvalue-space coefficient caching (dense: `msg.dot(V) * msg.dot(V_inv.T)`, sparse: per-site with multiplicity)
- [x] Analytical first and second derivatives (v1-only, v0 uses derivative-free Brent)
- [x] Newton's method in sqrt(t) space (`solve_newton`, safeguarded by a slope bracket and bisection, max 50 inner iterations, step in `[s - sqrt(s^2 + 1), s]`, absolute tolerance floor 1e-8, reduces indel Hessian singularity from O(1/t^2) to O(1/t))
- [x] Brent's method in sqrt(t) space (`solve_brent`, derivative-free via `argmin::BrentOpt`, default, matches v0 exactly)
- [x] Per-edge optimization method selection (`--opt-method`: `brent` default, `newton`)
- [x] Grid search for zero-boundary reconciliation (100 points, log-spaced grid with 0.5 subs/site minimum upper bound)
- [x] Zero branch length short-circuit (combined likelihood > 0.01 and derivative < 0 at zero)
- [x] `compute_derivatives` flag to skip derivative computation for log-likelihood-only evaluation
- [x] Typed aggregate likelihood results and convergence history (`LogLh`), with native `f64` derivative kernels
- [x] Collect dense contribution (`PartitionMarginalDense`)
- [x] Collect sparse contribution (`PartitionMarginalSparse`, multiplicity-weighted)
- [x] Unified mixed-partition evaluation (`evaluate_partitions()` sums metrics across partition types)
- [x] sqrt(t) reparameterization (Brent in s=sqrt(t) space, matching v0)
- [ ] Regularization penalty for profile-based optimization (v0: `exp(t^4/10000)` prevents unbounded growth)
- [ ] Hamming distance fallback when optimization fails (v0 falls back to observed distance)

## Convergence Loop

- [x] Iterative marginal reconstruction + optimization loop bounded by `--max-iter`
- [x] Early stop when absolute likelihood change is below `--dp`
- [x] Oscillation detection: stop when 2-step likelihood difference is below `--dp`
- [x] Worsened detection: stop and revert to best branch lengths when likelihood decreases from peak
- [x] Damping in marginal loop (`--damping`, `new*(1-d^i) + old*d^i` as in v0; default 0, v0 default 0.75: [kb/decisions/optimize-damping-disabled-by-default.md](../decisions/optimize-damping-disabled-by-default.md))
- [x] Damping floor (`DAMPING_FLOOR = 0.01`) prevents fully undamped late iterations when `--damping` is above 0
- [ ] Progressive per-iteration tolerance tightening (v0: `tol = 1e-8 + 0.01^(i+1)`, coarse early, tight late)
- [x] Bifurcating root special handling: v0 optimizes the combined length of the root children and keeps their ratio; v1 does the same for every chain of edges in series, single-child nodes included ([kb/decisions/optimize-series-edges-joint-sum.md](../decisions/optimize-series-edges-joint-sum.md))
- [ ] Convergence by sequence change count (v0 joint mode: stops when zero nucleotides change)

## Topology Cleanup (v1-only)

Once per iteration, `prune_and_merge_in_loop` simplifies the tree. Sparse partitions only; dense partitions carry no per-edge mutation lists, so the reversion and merge steps are inert under a dense-only run.

- [x] Zero-optimal edge collapse: internal edges the per-edge optimizer drove to exactly zero, carrying no substitutions or indels, are contracted into polytomies. Mirrors v0's `prune_short_branches` inside the loop.
- [x] Shared-mutation merge: siblings in a polytomy that carry identical substitutions are grouped under a new internal node (`merge_shared_mutation_branches`).
- [x] Reversion hoist: when a child edge reverts a substitution on the node's parent edge, a new node is inserted grouping that child with its sibling subtree, lifting the non-reverted substitutions above it. Removes one mutation per reverted position and adds none; branch lengths split proportionally to preserve root-to-node distances ([kb/decisions/optimize-polytomy-reversion-resolution.md](../decisions/optimize-polytomy-reversion-resolution.md)).
- [x] Helper-node retirement: mutation-free edges to nodes created during the routine are collapsed, dissolving the transient helpers the merge and hoist leave behind.
- [x] Per-step toggles: `--no-collapse-short-branches`, `--no-merge-siblings`, `--no-flip-parent-child` disable individual steps (all default on, so output is unchanged without them). Carried in-memory as `TopologyOps` on `OptimizeParams`. Merge and flip-parent-child are sparse-only, so both warn when disabled under `--dense` where they cannot run. When collapse is disabled the optimizer may leave zero-length internal edges in the output; they are kept as-is.

The merge, hoist, and retire steps form one per-polytomy routine (`resolve_polytomies`) run every iteration and driven to a fixpoint by a monotone (mutation count, node count) potential. A polytomy whose children all revert the same position collapses to a single reverting child that the two-children guard leaves in place ([kb/issues/M-optimize-reversion-hoist-single-child-residual.md](../issues/M-optimize-reversion-hoist-single-child-residual.md)).

## GTR Integration

- [x] `--model` and `--model-params` wired through `GtrModel` for all named models and inference
- [ ] GTR inference integrated into optimization loop (v0: `infer_gtr` parameter re-estimates model per iteration)
- [ ] GTR rate optimization (v0: `optimize_gtr_rate()` optimizes overall substitution rate mu)

## Reuse by Timetree

- [x] `PartitionOptimizeOps` trait implemented by both dense and sparse partitions
- [x] `collect_edge_contributions()` gathers contributions for one edge across partition types
- [x] `compute_branch_length_distribution()` evaluates log-likelihood (substitution + Poisson indel) on grid, converts to time-domain distribution
- [x] `evaluate_with_indels_log_lh_only()` for grid evaluation without derivatives (matches `optimize_edge_lengths()`'s edge evaluator)

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

- [ ] Automatic choice of the representation without `--dense`: `infer_dense()` always returns false, so the default is sparse
- [ ] Standalone `run_optimize_sparse()` zero-branch inconsistency
