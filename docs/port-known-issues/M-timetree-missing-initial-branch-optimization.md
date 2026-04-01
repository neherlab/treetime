# Timetree skips initial ML branch length optimization before time inference

v0's timetree calls `optimize_tree(max_iter=1)` before the first `make_time_tree()`, performing one round of ML branch-length optimization (Brent per-edge with damping). v1's timetree goes directly to distribution-based time inference (`run_timetree`) using input branch lengths from the tree file and clock regression, without an initial ML optimization pass. The ML optimization pre-step maximizes the phylogenetic likelihood <a id="cite-1"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[1](#ref-1)] before time inference begins, as described in the TreeTime pipeline <a id="cite-2"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[2](#ref-2)].

## v0 behavior

In `treetime.py`, before the iteration loop:

```python
# Line 243 (first call, with GTR inference)
self.optimize_tree(infer_gtr=infer_gtr, max_iter=1, method_anc=method_anc, **seq_kwargs)

# Line 266 (second call, after rerooting)
self.optimize_tree(max_iter=1, method_anc=method_anc, **seq_kwargs)
```

`optimize_tree(max_iter=1)` calls `optimize_tree_marginal(max_iter=1)` ([packages/legacy/treetime/treetime/treeanc.py#L1297-L1360](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1360)), which:

1. Runs marginal ancestral reconstruction
2. Executes one iteration of per-edge Brent optimization with damping (`damping=0.75`)
3. Re-runs ancestral reconstruction with updated branch lengths

This seeds `make_time_tree()` with ML-optimized branch lengths rather than raw input values.

Inside the iteration loop, v0 does NOT re-optimize branch lengths via `optimize_tree` in the normal path (only after polytomy resolution, and then with `max_iter=0` which skips the optimization loop).

## v1 behavior

`run_refinement_iteration()` ([packages/treetime/src/commands/timetree/refinement.rs#L21-L106](../../packages/treetime/src/commands/timetree/refinement.rs#L21-L106)) calls `update_marginal()` (ancestral reconstruction) and `run_timetree()` (distribution-based inference). Neither the pre-loop initialization nor the iteration loop calls `run_optimize_mixed()` or `initial_guess_mixed()`.

The infrastructure exists: `run_optimize_mixed()` is generic over `PartitionOptimizeOps`, and timetree partitions implement this trait. The code could be called but is not wired.

## Impact

Branch lengths entering the first time-tree pass are unoptimized (raw input values adjusted by clock regression). For well-resolved trees with accurate input branch lengths, the impact is small because the distribution-based inference finds reasonable node times regardless. For trees with poor initial branch lengths (e.g., from neighbor-joining or parsimony), the missing optimization pass may cause the time-tree inference to start from a worse initial state, affecting convergence.

## Proposed solution

Add one call to `run_optimize_mixed()` (or a subset of its logic) before the first `run_timetree()` in the timetree pipeline. This aligns with v0's `optimize_tree(max_iter=1)` and reuses the existing Newton/grid-search infrastructure.

Design decision: make the pre-step configurable (not hardcoded on/off). Trees from IQ-TREE or FastTree already have decent ML branch lengths and may not need a pre-step. Trees from neighbor-joining or parsimony benefit from it.

Considerations:

- The timetree command has its own iteration loop and convergence criteria. Adding ML optimization should not conflict with the distribution-based inference but needs testing for convergence behavior.
- v0's `optimize_tree_marginal` includes damping (`0.75`). The standalone `run_optimize_mixed()` does not have damping (no issue file yet for optimize oscillation/no-damping). A single undamped pass (`max_iter=1` equivalent) may be sufficient for initial seeding.
- The `initial_guess_mixed()` function computes branch lengths from mutation counts and could serve as a lighter alternative to full Newton optimization for the initial seeding step.
- A CLI flag (e.g. `--optimize-branch-lengths-pre-step` with auto/always/never modes, following the pattern of `--branch-length-initial-guess`) would allow users to control this.

## Implementation pointers

The insertion point is `run_timetree_estimation()` in `packages/treetime/src/commands/timetree/run.rs`, after partition initialization and before the first `run_timetree()` call. Timetree partitions already implement `PartitionOptimizeOps` (the trait required by `run_optimize_mixed`). The wiring is missing, not the infrastructure.

For the CLI flag, follow the `--branch-length-initial-guess` pattern: an enum with `auto`/`always`/`never` modes, defaulting to `auto`. In `auto` mode, run the pre-step when input branch lengths appear unoptimized (all equal, all zero, or high variance relative to substitution counts).

Consider whether `initial_guess_mixed()` is a sufficient lighter alternative when input branch lengths are entirely missing (no Newton needed, just substitution-count seeding).

## Cross-references

- [Command relationships](../reports/command-relationships/_index.md): Documents the ideal hierarchy where the optimize pre-step feeds timetree.
- [Iteration loop design](../reports/iterative-tree-refinement/9-iteration-loop.md): Loop structure and convergence criteria.

## References

1. <a id="ref-1"></a> Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-1)
2. <a id="ref-2"></a> Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-2)
