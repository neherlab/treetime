# Branch Length Optimization

v0 has no standalone optimize command. `TreeAnc.optimize_tree()` and `optimize_tree_marginal()` perform branch length optimization inline during timetree and ancestral workflows. v1 extracts this into a standalone `optimize` command and a reusable `PartitionOptimizeOps` trait system shared with timetree.

## Standalone Command (v1-only)

- [x] `optimize` command with `--tree`, `--alignment`, `--output-all`
- [x] Dense or sparse partition setup, one partition per `--alignment` file or per `--partition` entry, selected by `--representation` for the whole run
- [x] Several partitions share the branch lengths, with one relative rate per partition (`--partition-rates`, only with two or more partitions, default `estimate`: a rate step and a site-weighted normalization after each iteration; `equal` fixes every rate at 1; all rates normalized once before the first iteration when any partition has the model `infer`) and one indel rate per file in one combined Poisson term, with the indel events of a split file counted once ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md))
- [x] Short-branch floors use the total length of all files
- [x] Tips without a sequence in a file are filled as unknown; the run stops when a file misses more than a third of the tips
- [x] Initial branch-length guess from observed mutation counts (excludes deletion and ambiguous positions)
- [x] `--initial-branch-lengths` flag: `auto` (selective fill, treats zero BL as invalid when indels present), `always` (overwrite all), `never` (error on missing/NaN, and on zero branch length when indels are present)
- [x] Output Newick and Nexus trees with the branch mutations (Nexus always, Newick with `--newick-annotations beast` or `nhx`)
- [x] Model document (`--output-model`): the alignment files with their number, path, sites and indel rate, and one model entry per partition with its number, name, alignment number, sites, state characters
- [x] Augur-compatible node data JSON (`optimize.node-data.json`, `--output-node-data`): per-node optimized `branch_length` plus `alignment`/`input_tree` metadata, equivalent to `augur refine` without `--timetree`

## Rerooting (v1-only, date-free)

The optimize command has no sampling dates, so only date-free rooting policies are available. Rerooting is opt-in; the input root is kept by default.

- [x] `--reroot=min-dev`: minimize the variance of root-to-tip distances with the edge search at a fixed rate of zero, every tip counted ([kb/algo/reroot.md](../algo/reroot.md)). This implements v0's documented fixed-zero-slope meaning; actual v0 candidate scoring retains an estimated-rate term ([kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md)).
- [x] `--reroot=midpoint`: root at the middle of the longest tip-to-tip path, as `augur refine --root mid_point` (Biopython `Tree.root_at_midpoint()`) does ([kb/algo/reroot.md](../algo/reroot.md), "Midpoint Rooting")
- [x] `--outgroup=A,B,...`: root at the middle of the branch above the MRCA of the named tips in the unrooted tree, no scoring; a repeated reroot keeps the root ([kb/decisions/reroot-named-node-at-unrooted-branch-middle.md](../decisions/reroot-named-node-at-unrooted-branch-middle.md)).
- [x] `--keep-root`: explicit form of the default (keep the input root); mutually exclusive with the reroot options.
- [x] Date-dependent methods (`least-squares`, `oldest`, `clock-filter`) rejected at argument validation
- [x] Reroot once right after the initial guess, on the starting branch lengths, then partition reconciliation and branch-length optimization in the main loop ([kb/decisions/reroot-date-free-methods-root-once.md](../decisions/reroot-date-free-methods-root-once.md))

## Per-Edge Likelihood (v0 parity via different method)

v0 uses Brent's method (`scipy.optimize.minimize_scalar`) in sqrt(t) space with Hamming distance bracket. v1 offers two methods via `--optimizer`: Brent and Newton-Raphson, both in $\sqrt{t}$ space. Default is `brent` (matches v0). Both v0 and v1 use eigendecomposition-based likelihood (`expQt = V diag(exp(lambda*t)) V_inv`).

- [x] Eigenvalue-space coefficient caching (dense: `msg.dot(V) * msg.dot(V_inv.T)`, sparse: per-site with multiplicity)
- [x] Analytical first and second derivatives (v1-only, v0 uses derivative-free Brent)
- [x] Newton's method in sqrt(t) space (`solve_newton`, safeguarded by a slope bracket and bisection, max 50 inner iterations, step in `[s - sqrt(s^2 + 1), s]`, absolute tolerance floor 1e-8, reduces indel Hessian singularity from O(1/t^2) to O(1/t))
- [x] Brent's method in sqrt(t) space (`solve_brent`, derivative-free via `argmin::BrentOpt`, default, matches v0 exactly)
- [x] Per-edge optimization method selection (`--optimizer`: `brent` default, `newton`)
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
- [x] Early stop when absolute likelihood change is below `--tolerance`
- [x] Oscillation detection: stop when 2-step likelihood difference is below `--tolerance`
- [x] Worsened detection: stop and revert to best branch lengths when likelihood decreases from peak
- [x] Damping in marginal loop (`--damping`, `new*(1-d^i) + old*d^i` as in v0; default 0, v0 default 0.75: [kb/decisions/optimize-damping-disabled-by-default.md](../decisions/optimize-damping-disabled-by-default.md))
- [x] Damping floor (`DAMPING_FLOOR = 0.01`) prevents fully undamped late iterations when `--damping` is above 0
- [ ] Progressive per-iteration tolerance tightening (v0: `tol = 1e-8 + 0.01^(i+1)`, coarse early, tight late)
- [x] Bifurcating root special handling: v0 optimizes the combined length of the root children and keeps their ratio; v1 does the same for every chain of edges in series, single-child nodes included ([kb/decisions/optimize-series-edges-joint-sum.md](../decisions/optimize-series-edges-joint-sum.md))
- [ ] Convergence by sequence change count (v0 joint mode: stops when zero nucleotides change)

## Topology Cleanup (v1-only)

Once per iteration, `fn prune_and_merge_in_loop()` simplifies the tree, and once after the loop `fn collapse_short_internal_edges()` removes the short internal edges. The merge, hoist, and retirement steps read per-edge mutation lists, which only sparse partitions carry, so they are inert under a dense-only run.

- [x] Zero-optimal edge collapse: every internal edge that the per-edge optimizer drove to exactly zero is contracted into a polytomy, also when it carries Fitch substitutions or indels, because the transition matrix at length zero is the identity and the removal leaves the substitution likelihood unchanged. An edge stays only when its indel term enters the likelihood and needs a positive length (not under `--no-indels`). The Fitch substitutions, marginal estimates, and indels of the edge move to its child edges; a composed Fitch substitution at a site where the child subtree has no data is dropped
- [x] Short-branch collapse after the loop: on the final branch lengths, every internal edge of length zero, and every internal edge shorter than $0.1/L$ ($L$ is the number of sites of all partitions) whose marginal estimate carries no substitution in any partition, is contracted; no merge or hoist follows. An edge whose indel term needs a positive length stays. This is v0's marginal-mode `prune_short_branches` after `optimize_tree`: its test $P(0) > 0.1$ holds exactly when the parent and child sequences agree at every compared site
- [x] Shared-mutation merge: siblings in a polytomy that share substitutions are grouped under a new internal node (`fn merge_single_polytomy()`). A substitution counts only when Fitch and the marginal reconstruction both place it on the edge (`fn SparseReconstruction::supported_fitch_subs()`); an edge without a marginal estimate, as in `prune`, keeps all its Fitch substitutions. Indels count only when they enter the likelihood, so `--no-indels` merges on substitutions only. A group of all children of a node is never merged, because it would leave a single-child node
- [x] Merge ledger: a merge stays only while the optimizer keeps its new edge above length zero. When the collapse removes the edge into a helper node of a merge, `struct MergeLedger` records the set of its children as rejected, and no later merge of the run groups the same children, under any parent. `prune` starts with an empty ledger
- [x] Reversion hoist: when a child edge reverts a substitution on the node's parent edge, a new node is inserted grouping that child with its sibling subtree, lifting the non-reverted substitutions above it. Removes one mutation per reverted position and adds none; branch lengths split proportionally to preserve root-to-node distances. A reversion counts only when Fitch and the marginal reconstruction share both substitutions. When the hoist leaves the node with one child, its parent edge is collapsed, which cancels the last reversion ([kb/decisions/optimize-polytomy-reversion-resolution.md](../decisions/optimize-polytomy-reversion-resolution.md))
- [x] Helper-node retirement: mutation-free edges to nodes created during the routine are collapsed, dissolving the transient helpers the merge and hoist leave behind. Under `--no-indels` an edge with only indels is mutation-free
- [x] Sparse data of new nodes: a node that a merge or a hoist creates takes its gap, unknown, and non-character ranges from its children, and its composition from the sequence of the node above it with the substitutions of its new edge (`fn new_node_obs()`)
- [x] Per-step toggles: `--no-collapse-short-branches`, `--no-merge-shared-mutations`, `--no-hoist-reversions` disable individual steps (all default on, so output is unchanged without them). Carried in-memory as `TopologyOps` on `OptimizeParams`. `--no-collapse-short-branches` disables both collapses. Merge and flip-parent-child are sparse-only, so `--no-merge-shared-mutations` and `--no-hoist-reversions` together with `--representation dense` stop the run with an argument error. When collapse is disabled the optimizer may leave zero-length internal edges in the output; they are kept as-is

The merge, hoist, and retire steps form one per-polytomy routine (`fn resolve_polytomies()`) run every iteration and driven to a fixpoint by a monotone (mutation count, node count) potential.

- [x] Stop rule on one topology: `struct TopologyHistory` clears the log-likelihood history and the best iteration on every topology change, so the convergence, oscillation, and decrease checks compare iterations of the current topology only. A non-finite log-likelihood stops the loop at once
- [x] Best iteration on one topology: the loop returns the branch lengths, partition rates, and indel rates of the best iteration of the final topology, and never records the iteration directly after a topology change, whose new edges still have their start lengths
- [x] Final branch-length step: when the iteration limit ends the loop directly after a topology change, one more branch-length optimization and marginal update run without topology cleanup. Every edge of the result is then optimized on the final topology at least once, except that the post-loop collapse adds lengths below $0.1/L$ to child edges

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

## Short Branch Handling

- [x] Zero-branch-length derivative check before optimization
- [x] Short branch pruning after optimization, as v0 marginal mode after its loop (`bl < 0.1 * one_mutation` and P(zero) > 0.1); see Topology Cleanup
- [ ] MIN_BRANCH_LENGTH floor for GTR calculations (v0: `1e-3 * one_mutation`)

## Current Limitations

- [ ] Automatic choice of the representation with `--representation auto` (the default): `infer_dense()` always returns false, so the default is sparse
- [ ] Standalone `run_optimize_sparse()` zero-branch inconsistency
