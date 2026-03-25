# Chapter 10: Implementation status and recommendations

[Back to index](_index.md) | Previous: [Chapter 9: The iteration loop](9-iteration-loop.md)

## Feature matrix

### Topology cleanup operations

| Operation                      | v0 optimize      | v0 timetree       | v1 optimize     | v1 prune                 | v1 timetree     | RAxML     | IQ-TREE   |
| ------------------------------ | ---------------- | ----------------- | --------------- | ------------------------ | --------------- | --------- | --------- |
| Zero-length detection          | Threshold + prob | Via optimize_tree | Derivative sign | N/A                      | N/A             | Min clamp | Min clamp |
| Zero-length pruning            | Inside loop      | After resolution  | Not implemented | --prune-short            | Not implemented | N/A       | N/A       |
| Greedy polytomy resolution     | N/A              | Inside loop       | Not implemented | Not implemented          | Inside loop     | N/A       | N/A       |
| Stochastic polytomy resolution | N/A              | Inside loop       | N/A             | N/A                      | Not implemented | N/A       | N/A       |
| Shared-mutation merging        | N/A              | N/A               | Not implemented | --merge-shared-mutations | N/A             | N/A       | N/A       |

### Optimization methods

| Aspect                      | v0                          | v1                            | RAxML-NG                     | IQ-TREE                        |
| --------------------------- | --------------------------- | ----------------------------- | ---------------------------- | ------------------------------ |
| Per-edge optimizer          | Brent in sqrt(t)            | Newton-Raphson + grid         | Newton-Raphson + clamp       | NR + bisection                 |
| Outer-loop damping          | 0.75 exponential (marginal) | 0.75 exponential (all)        | None (per-branch safeguards) | None (stochastic perturbation) |
| Post-topology branch re-opt | Yes                         | No                            | N/A                          | N/A                            |
| Convergence                 | N_diff or delta_LH          | delta_LH or n_diff+n_resolved | LH < epsilon                 | LH plateau                     |

### CLI flags

| Flag                        | v0 optimize          | v0 timetree            | v1 optimize        | v1 prune        | v1 timetree                                     |
| --------------------------- | -------------------- | ---------------------- | ------------------ | --------------- | ----------------------------------------------- |
| prune_short / --prune-short | True (default)       | False after resolution | Not implemented    | Yes (threshold) | Not implemented                                 |
| --merge-shared-mutations    | N/A                  | N/A                    | Not implemented    | Yes             | N/A                                             |
| --resolve-polytomies        | N/A                  | True (default)         | N/A                | N/A             | Yes                                             |
| --keep-polytomies           | N/A                  | Available              | N/A                | N/A             | Parsed, not wired (`N-timetree-dead-cli-flags`) |
| --stochastic-resolve        | N/A                  | Available              | N/A                | N/A             | Not implemented                                 |
| --damping                   | N/A (hardcoded 0.75) | N/A                    | Yes (default 0.75) | N/A             | N/A                                             |

## Defects by severity

### Blocking

**Edge collapse uses substitution union instead of composition.** `collapse_sparse_edge()` merges substitutions via `iterator_union()`. The Markov semigroup property requires composition ([Chapter 2](2-substitution-models.md)): `A->G` then `G->T` = net `A->T`, not `{A->G, G->T}`. The existing `compose_substitutions()` function is the correct replacement. Tracked: `M-prune-collapse-uses-union-not-composition`. Blocks all topology cleanup integration.

### Important

**No topology cleanup in optimize loop.** v0 prunes zero-length branches inside each iteration. v1 detects zero-optimal branches but does not collapse them. Zero-length edges accumulate, wasting computation and preventing polytomy resolution. Tracked: `M-optimize-no-topology-cleanup-in-loop`. See [Chapter 9](9-iteration-loop.md) for the proposed loop.

**No branch re-optimization after timetree polytomy resolution.** v0 calls `optimize_tree(max_iter=0)` after resolving polytomies. v1 skips this. The time inference adjusts branch lengths indirectly through distributions, which partially compensates.

### Minor

**`merge_compressed` parameter unused.** The stretched/compressed distinction from v0 is declared but not implemented in v1. All children are merged regardless.

**Shared-mutation merging ignores shared indels.** Indels on child edges are preserved without merging identical indels to the parent edge. Conservative behavior.

**O(n^2) complexity for large polytomies.** Both greedy algorithms use pairwise comparison. Impractical for SARS-CoV-2-scale polytomies.

**Stochastic polytomy resolution not implemented.** Tracked: `N-timetree-stochastic-polytomy-unimplemented`. Greedy mode is deprecated in v0.

## Positive points

**Zero-length detection is mathematically sound.** v1's derivative-sign method is the local optimality condition at the boundary (exact for unimodal models such as JC69). See [Chapter 6](6-zero-length-branches.md). No heuristic thresholds.

**Damping matches v0 and SOR theory.** Exponential damping with `d=0.75` (Young 1954). See [Chapter 9](9-iteration-loop.md).

**Shared-mutation merging is correct and well-tested.** 515 lines of tests, smoke tests. Matches the design document specification.

**Timetree polytomy resolution uses proper likelihood criterion.** Brent optimization with `zero_branch_slope = mu * L`, matching v0.

**Convergence monitoring is more sophisticated.** `TimetreeOptimizer` tracks multiple metrics.

**Partition reconciliation after topology changes.** `reconcile_topology()` ensures data consistency for new nodes and edges.

## Recommendations

### R1 (blocking): Fix composition bug

Replace `iterator_union()` with `compose_substitutions()` in `collapse_sparse_edge()`. One function call change. Prerequisite for everything else.

### R2 (high): Integrate topology cleanup into optimize loop

Implement `prune_short_branches()` with v0's compound criterion. Integrate after damping, before next E-step. See [Chapter 9](9-iteration-loop.md) for the loop structure. Add `--prune-short` flag to optimize command (default enabled).

### R3 (high): Extract topology operations to shared module

Move `merge_shared_mutation_branches()` and `collapse_sparse_edge()` to a shared module callable from both optimize and prune commands.

### R4 (medium): Add `--merge-shared-mutations` to optimize command

Default disabled (new feature not in v0).

### R5 (low): Implement stochastic polytomy resolution

For SARS-CoV-2-scale datasets. Priority depends on target dataset scale.

## Ledger cross-references

| Issue                     | File                                           | Status |
| ------------------------- | ---------------------------------------------- | ------ |
| Composition bug           | `M-prune-collapse-uses-union-not-composition`  | Open   |
| No topology cleanup       | `M-optimize-no-topology-cleanup-in-loop`       | Open   |
| Stochastic resolution     | `N-timetree-stochastic-polytomy-unimplemented` | Open   |
| Polytomy numerical issues | `N-timetree-polytomy-numerical-robustness`     | Open   |
| Optimize loop extraction  | `L-optimize-loop-not-extracted`                | Open   |

## References

- Young, D. M. 1954. "Iterative Methods for Solving Partial Difference Equations." _Trans. Amer. Math. Soc._ 76(1):92-111. https://doi.org/10.1090/S0002-9947-1954-0059635-7
- Steel, M., and D. Penny. 2000. "Parsimony, Likelihood, and the Role of Models." _Mol. Biol. Evol._ 17(6):839-850. https://doi.org/10.1093/oxfordjournals.molbev.a026364
- Sagulenko, P., V. Puller, and R. A. Neher. 2018. "TreeTime." _Virus Evol._ 4(1):vex042. https://doi.org/10.1093/ve/vex042
