# Zero branch length optimization

How TreeTime v1 decides when an edge has zero optimal branch length, handles boundary conditions, and collapses degenerate edges. The central files are [optimize/dispatch.rs](../../packages/treetime/src/optimize/dispatch.rs) and [optimize/zero_boundary.rs](../../packages/treetime/src/optimize/zero_boundary.rs).

## Model classification

Whether the zero-branch shortcut applies depends on the substitution model. `fn GTR::certified_decay_rate()` ([gtr/gtr.rs#L121-L123](../../packages/treetime/src/gtr/gtr.rs#L121-L123)) records this classification: it returns the decay rate shared by every nonzero eigenvalue, and `None` for other models. `fn GTR::new()` sets the unit decay rate when all exchangeabilities are equal and positive ([gtr/gtr.rs#L100-L101](../../packages/treetime/src/gtr/gtr.rs#L100-L101), `fn has_equal_positive_exchangeabilities()` [gtr/gtr.rs#L222](../../packages/treetime/src/gtr/gtr.rs#L222)).

JC69 (`fn jc69()` [gtr/get_gtr.rs#L120](../../packages/treetime/src/gtr/get_gtr.rs#L120)), F81 (`fn f81()` [gtr/get_gtr.rs#L169](../../packages/treetime/src/gtr/get_gtr.rs#L169)), and binary models (equal exchangeabilities by construction) are classified as <a id="gloss-use-1"></a>unimodal <sup>[1](#gloss-1)</sup>. K80, HKY85, TN93, and general GTR are not. The theoretical basis is <a id="cite-1a"></a>[Dinh and Matsen 2017](https://doi.org/10.1214/16-AAP1240) [[1](#ref-1)], Corollary 3.1: models with a single distinct nonzero eigenvalue have at most one stationary point on $(0, \infty)$.

This classification drives the entire zero-branch pipeline. Unimodal models get a cheap derivative test. Non-unimodal models skip the test and rely on post-solver grid search to catch boundary errors.

## Per-edge optimization pipeline

Each edge goes through three phases inside `fn optimize_edge()` ([optimize/dispatch.rs#L272-L308](../../packages/treetime/src/optimize/dispatch.rs#L272-L308), called by `run_optimize_mixed_inner()` [optimize/dispatch.rs#L48](../../packages/treetime/src/optimize/dispatch.rs#L48)): input protection, the optimizer itself, and post-optimization reconciliation.

### Phase 1: Input domain protection

Before the optimizer runs, two guards ensure the starting branch length is in the evaluator's well-defined domain.

**Indel-bearing edges.** The Poisson derivative diverges at $t = 0$ when $k > 0$, where $t$ is the branch length and $k$ is the indel count on the edge. The starting point is bumped to $\max(k / \text{rate}, 1/L)$ (or $1/L$ when the indel rate is zero), where $L$ is the total alignment length across all partitions ([optimize/dispatch.rs#L280-L286](../../packages/treetime/src/optimize/dispatch.rs#L280-L286)).

**Invalid-at-zero edges.** When any site has $L_i(0) \leq 0$ (where $L_i$ is the per-site likelihood, computed as $L_i(0) = \sum_c k_{ic}$ with $k_{ic}$ being the eigenvalue-space coefficients), evaluating $\ln L_i$ produces $-\infty$. The starting point is bumped to $1/L$ ([optimize/dispatch.rs#L288-L290](../../packages/treetime/src/optimize/dispatch.rs#L288-L290)).

During optimization, `min_branch_length_for_indels()` ([optimize/zero_boundary.rs#L13-L19](../../packages/treetime/src/optimize/zero_boundary.rs#L13-L19)) provides a floor for Newton/Brent steps: $0.01/L$ for indel-bearing edges, $0$ otherwise.

### Phase 2: Pre-dispatch shortcut

Before dispatching to a solver, `is_zero_branch_optimal()` ([optimize/zero_boundary.rs#L21-L40](../../packages/treetime/src/optimize/zero_boundary.rs#L21-L40)) checks whether the edge can be set to zero immediately ([optimize/dispatch.rs#L292-L294](../../packages/treetime/src/optimize/dispatch.rs#L292-L294)). This requires four conditions to hold simultaneously:

1. Every partition uses a unimodal model with the same certified decay rate (`certified_decay_rate()`, [partition/optimize/contribution.rs#L72-L74](../../packages/treetime/src/partition/optimize/contribution.rs#L72-L74); `shared_certified_decay_rate()`, [optimize/zero_boundary.rs#L120-L127](../../packages/treetime/src/optimize/zero_boundary.rs#L120-L127))
2. No indels on this edge (`indel_count == 0`)
3. Every site has positive finite likelihood at $t = 0$ (`all_sites_valid_at_zero()`, [partition/optimize/contribution.rs#L58-L63](../../packages/treetime/src/partition/optimize/contribution.rs#L58-L63))
4. The total derivative of log-likelihood at $t = 0$ is finite and negative

The derivative is computed by `zero_branch_length_derivative()` ([partition/optimize/contribution.rs#L76-L86](../../packages/treetime/src/partition/optimize/contribution.rs#L76-L86)):

$$\frac{d}{dt} \log L(0) = \sum_i m_i \frac{\sum_c k_{ic} \lambda_c}{\sum_c k_{ic}}$$

where $m_i$ is the multiplicity of site pattern $i$ and $\lambda_c$ are the GTR eigenvalues. This is a weighted average of eigenvalues per site, summed over all sites.

If the derivative is negative, the log-likelihood decreases as $t$ moves away from zero, and unimodality guarantees zero is the global maximum. The edge is set to zero and the solver is skipped. If any condition fails, the edge falls through to Newton/Brent.

### Phase 3: Post-dispatch reconciliation

After the solver returns a candidate branch length, `reconcile_zero_boundary()` ([optimize/zero_boundary.rs#L42-L65](../../packages/treetime/src/optimize/zero_boundary.rs#L42-L65)) checks for three failure modes:

1. Positive but worse than zero. The solver returned a tiny positive value (e.g. $10^{-12}$) when zero has higher likelihood. Brent's method places a strictly positive lower bound on the bracket, so it cannot return exact zero even when zero is optimal. Detected by `is_zero_better_than_grid_best()` ([optimize/zero_boundary.rs#L107-L118](../../packages/treetime/src/optimize/zero_boundary.rs#L107-L118)), which compares the log-likelihood at zero against $\log L(t_{\text{best}})$ at the best positive candidate.

2. Exactly zero on a non-unimodal model. Newton/Brent clamped to zero via step limiting, but a better positive mode may exist. The K80 counterexample from <a id="cite-1b"></a>[Dinh and Matsen 2017](https://doi.org/10.1214/16-AAP1240) [[1](#ref-1)] with transition/transversion ratio $\kappa = 3$ has a local maximum at $t \approx 0.167$, a local minimum at $t \approx 1.25$, and its global maximum at $t \approx 6.37$ (values from the exact coefficients of the test fixture, `fn make_dinh_matsen_k80_contribution()` in `optimize/__tests__/test_dispatch_zero_boundary.rs`). A solver that clamps to zero from a starting point $t_0 \approx 0.6$ misses the local maximum. The reconcile step finds the local maximum only because it lies inside the grid; the global maximum lies outside the grid and outside v0's search range $t \leq 4$, and the tests of this path assert a positive point that beats zero within the grid, not the global maximum.

3. Exactly zero with invalid evaluation domain. The solver returned zero on a partition where $L_i(0) \leq 0$. The caller bumped the input away from zero (Phase 1), but step clamping brought it back. Preserving this zero would re-enter the invalid domain on the next marginal update.

All three cases route to `grid_search_inner()` ([optimize/zero_boundary.rs#L67-L95](../../packages/treetime/src/optimize/zero_boundary.rs#L67-L95)), which evaluates 100 log-spaced candidates produced by `grid_search_branch_lengths()` ([optimize/zero_boundary.rs#L97-L105](../../packages/treetime/src/optimize/zero_boundary.rs#L97-L105)). The grid spans from $0.1/L$ to $\max(1.5t + 1/L, 0.5)$, with log spacing providing uniform resolution per decade across the 3-4 orders of magnitude that branch lengths span. The grid search compares its best positive candidate against zero and returns whichever has higher log-likelihood.

## Topology cleanup

After branch-length optimization, the optimize loop collapses zero-length edges and restructures the tree.

`find_zero_optimal_internal_edges()` ([optimize/run_loop.rs#L299](../../packages/treetime/src/optimize/run_loop.rs#L299)) runs after `run_optimize_mixed_inner()` and before damping ([optimize/run_loop.rs#L230](../../packages/treetime/src/optimize/run_loop.rs#L230)). It collects internal edges with branch length exactly zero and no mutations or indels on any partition. Leaf edges are excluded, matching v0's `prune_short_branches` which skips terminals.

`prune_and_merge_in_loop()` ([optimize/run_loop.rs#L325](../../packages/treetime/src/optimize/run_loop.rs#L325)) then overrides damped branch lengths back to zero for these edges and collapses each one ([optimize/run_loop.rs#L359](../../packages/treetime/src/optimize/run_loop.rs#L359)). The underlying `collapse_edge()` ([optimize/topology/collapse.rs#L7-L43](../../packages/treetime/src/optimize/topology/collapse.rs#L7-L43)) removes the target node, reparents its children to the source, sums branch lengths, composes substitutions via the Markov semigroup property (reversions cancel), concatenates indels, and drops stale partition entries. `collapse_edge()` is also used by the prune command ([prune/prune.rs#L64](../../packages/treetime/src/prune/prune.rs#L64), [prune/prune.rs#L121](../../packages/treetime/src/prune/prune.rs#L121)).

After collapsing, sibling branches in newly formed polytomies that share identical substitutions are merged under new internal nodes.

## Data flow

```mermaid
graph TD
    subgraph pre ["Pre-dispatch (unimodal models only)"]
        valid["all_sites_valid_at_zero()"]
        unimodal["certified_decay_rate()"]
        deriv["zero_branch_length_derivative()"]
        is_zero["is_zero_branch_optimal()"]
        valid --> is_zero
        unimodal --> is_zero
        deriv --> is_zero
    end

    subgraph solver ["Solver dispatch"]
        newton["Newton"]
        brent["Brent"]
    end

    subgraph post ["Post-dispatch reconciliation"]
        reconcile["reconcile_zero_boundary()"]
        grid["grid_search_inner()"]
        zero_better["is_zero_better_than_grid_best()"]
        reconcile --> grid
        grid --> zero_better
    end

    subgraph collapse ["Topology cleanup (optimize loop)"]
        find["find_zero_optimal_internal_edges()"]
        prune["prune_and_merge_in_loop()"]
        coll["collapse_edge()"]
        find --> prune
        prune --> coll
    end

    run["run_optimize_mixed_inner()"]
    run -->|"indel_count == 0"| is_zero
    is_zero -->|"true: set bl=0, skip solver"| done["edge.set_branch_length(0)"]
    is_zero -->|"false"| solver
    solver --> reconcile
    reconcile -->|"final branch length"| done

    loop["run_optimize_loop()"]
    loop --> run
    run --> find
```

## Decision matrix

| Model class                          | Indels  | Sites valid at zero | Path                                                                                          |
| ------------------------------------ | ------- | ------------------- | --------------------------------------------------------------------------------------------- |
| Unimodal (JC69, F81, binary)         | none    | yes                 | Derivative shortcut. Negative derivative -> zero. Otherwise solver.                           |
| Unimodal                             | none    | no                  | Solver with bumped start ($1/L$). Post-dispatch reconciliation if solver returns zero.        |
| Unimodal                             | present | any                 | Solver with Poisson floor. Zero never optimal (Poisson $-\infty$).                            |
| Non-unimodal (K80, HKY85, TN93, GTR) | none    | yes                 | Solver. Post-dispatch grid search if solver returns zero or a positive value worse than zero. |
| Non-unimodal                         | none    | no                  | Solver with bumped start. Post-dispatch grid search if solver returns zero (invalid domain).  |
| Non-unimodal                         | present | any                 | Solver with Poisson floor. Zero never optimal.                                                |

## Glossary

1. <a id="gloss-1"></a> **Unimodal (branch-length likelihood).** A per-edge log-likelihood function with at most one local maximum on $(0, \infty)$. Proven for models with a single distinct nonzero eigenvalue: JC69, F81, and binary symmetric models ([Dinh and Matsen 2017](https://doi.org/10.1214/16-AAP1240) [[1](#ref-1)], Corollary 3.1). Models with multiple distinct nonzero eigenvalues (K80, HKY85, TN93, GTR) can have two local maxima. [↩](#gloss-use-1)

## References

1. <a id="ref-1"></a> Dinh, Vu, and Frederick A. Matsen IV. 2017. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." _Annals of Applied Probability_ 27(3):1646-1677. https://doi.org/10.1214/16-AAP1240 [↩¹](#cite-1a) [↩²](#cite-1b)
