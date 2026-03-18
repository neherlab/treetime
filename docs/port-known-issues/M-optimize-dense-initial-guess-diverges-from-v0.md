# Dense initial branch length guess diverges from v0

The dense partition's initial branch length estimate in `initial_guess_mixed()` uses a hard argmax Hamming distance on per-edge messages. v0 uses a soft Hamming distance (profile dot product) on full marginal profiles. Neither the current v1 criterion nor the pre-PR #470 criterion matches v0.

## Problem

`initial_guess_mixed()` ([packages/treetime/src/commands/optimize/optimize_unified.rs#L252-L286](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L252-L286)) computes per-edge branch lengths as the fraction of positions with detected substitutions. For dense partitions, substitution detection is delegated to `PartitionMarginalDense::edge_subs()` ([packages/treetime/src/representation/partition/marginal_dense.rs#L67-L96](../../packages/treetime/src/representation/partition/marginal_dense.rs#L67-L96)), which compares `argmax_first` of per-edge message profiles:

```rust
let parent_state = self.alphabet.char(argmax_first(&parent).unwrap_or(0));
let child_state = self.alphabet.char(argmax_first(&child).unwrap_or(0));
if parent_state != child_state {
  subs.push(Sub::new(parent_state, pos, child_state)?);
}
```

This is a hard binary comparison on per-edge messages (`edge.msg_to_parent.dis`, `edge.msg_to_child.dis`). Each position contributes either 0 or 1 to the difference count.

## v0 reference

v0's marginal branch length optimization in `optimize_tree_marginal()` ([packages/legacy/treetime/treetime/treeanc.py#L1297-L1346](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1346)) calls `optimal_marginal_branch_length()` which calls `marginal_branch_profile(node)` to obtain full marginal profiles, then passes them to `optimal_t_compressed(profiles=True)`.

### v0 profiles

`marginal_branch_profile(node)` ([packages/legacy/treetime/treetime/treeanc.py#L1122-L1146](../../packages/legacy/treetime/treetime/treeanc.py#L1122-L1146)) returns two profiles at the child node:

- `pp = node.marginal_outgroup_LH` - likelihood contribution from everything above this node (parent side), propagated down to the child
- `pc = node.marginal_subtree_LH` - likelihood contribution from the subtree below this node (child side)

These are full marginal likelihood vectors (shape `[L, n_states]`), computed after the backward and forward passes. The full posterior at the child is proportional to `pp * pc`. The two components represent complementary information about the child's state: what the rest of the tree says vs what the child's own subtree says.

### v0 soft Hamming distance

`optimal_t_compressed(profiles=True)` ([packages/legacy/treetime/treetime/gtr.py#L871-L876](../../packages/legacy/treetime/treetime/gtr.py#L871-L876)) computes:

```python
hamming_distance = 1 - np.sum(multiplicity * np.sum(pp * pc, axis=1)) / np.sum(multiplicity)
```

This is a weighted average of `1 - dot(pp[i], pc[i])` across alignment positions. For each position, `dot(pp[i], pc[i])` measures the overlap between the parent-side and child-side profiles:

- Sharp profiles (one state dominates): `dot ≈ 1` when both agree on the state, `dot ≈ 0` when they disagree. Approximates hard Hamming.
- Uncertain profiles (flat/near-uniform): for a 4-state nucleotide alphabet, uniform profiles give `dot = sum(0.25 * 0.25 * 4) = 0.25`, so each uncertain position contributes `1 - 0.25 = 0.75` to the distance. This reflects the information that the branch is long enough to erase the phylogenetic signal at this position.

### v0 role of Hamming distance

The soft Hamming distance serves as the bracket midpoint for Brent's method ([packages/legacy/treetime/treetime/gtr.py#L881-L891](../../packages/legacy/treetime/treetime/gtr.py#L881-L891)):

```python
opt = minimize_scalar(
    _neg_prob,
    bracket=[-sqrt(MAX_BRANCH_LENGTH), sqrt(hamming_distance), sqrt(MAX_BRANCH_LENGTH)],
    args=(seq_pair, multiplicity),
    method='brent',
)
new_len = opt['x'] ** 2
```

v0 parameterizes branch length as `t^2` (enforcing non-negativity). The bracket midpoint `sqrt(hamming_distance)` starts Brent near the expected optimum. Brent is a bracketed method: it always finds the minimum within the bracket regardless of the midpoint, but a good midpoint accelerates convergence.

### v0 full optimization objective

Brent minimizes `-prob_t_profiles` ([packages/legacy/treetime/treetime/gtr.py#L922-L963](../../packages/legacy/treetime/treetime/gtr.py#L922-L963)), which computes the log-likelihood of observing the parent-child profile pair at distance `t`:

```python
Qt = self.expQt(t)
res = np.einsum('ai,ij,aj->a', pc, Qt, pp)
logP = np.sum(multiplicity * np.log(res + SUPERTINY))
```

This marginalizes over all possible parent and child states, weighted by the GTR transition matrix `P(t) = exp(Q*t)`. The Hamming distance is just the starting point; the full model-based likelihood determines the final branch length.

## v1 history

Three versions of the dense initial guess have existed in v1:

### v1-original (before PR #470)

Inline code in `initial_guess_mixed()`:

```rust
for (prof_parent, prof_child) in izip!(rows_parent, rows_child) {
  if prof_parent[prof_child.argmax().unwrap()] < 0.5 {
    differences += 1;
  }
}
```

This checked whether the parent's probability at the child's most likely state was below 0.5. For a 4-state alphabet:

- Uncertain positions (near-uniform, each state ≈ 0.25): parent probability at child's argmax ≈ 0.25 < 0.5. Always counts as a difference. Behavior: overcounts relative to hard Hamming.
- Sharp agreeing profiles: parent probability at shared state is high (> 0.5). No difference counted.
- Sharp disagreeing profiles: parent probability at child's best state is low (< 0.5). Difference counted.

### v1-PR470 (PR #470, commit `7d6a14f8`)

Extracted to `PartitionMarginalDense::edge_subs()` as part of `PartitionOptimizeOps` trait unification. Changed to argmax comparison:

```rust
let parent_state = self.alphabet.char(argmax_first(&parent).unwrap_or(0));
let child_state = self.alphabet.char(argmax_first(&child).unwrap_or(0));
if parent_state != child_state { ... }
```

This is hard Hamming on per-edge message argmax states. For a 4-state alphabet:

- Uncertain positions: argmax is deterministic but arbitrary (`argmax_first` picks the first maximum). Two independent near-uniform profiles may or may not have the same argmax depending on noise. Behavior: non-deterministic difference count for uncertain positions.
- Sharp profiles: same as v0 hard Hamming.

### v1-current (after PR #470, commit `65ad978e`)

Added gap filtering via `edge_effective_length()` and graph parameter for gap range lookup. The counting criterion remained argmax comparison.

## Divergence analysis

The three approaches (v0 soft, v1-original threshold, v1-current argmax) differ at uncertain positions. For sharp profiles, all three give equivalent results.

| Scenario                            | v0 soft Hamming | v1-original threshold | v1-current argmax     |
| ----------------------------------- | --------------- | --------------------- | --------------------- |
| Sharp, agree                        | ≈ 0             | 0                     | 0                     |
| Sharp, disagree                     | ≈ 1             | 1                     | 1                     |
| Uniform (4-state)                   | 0.75            | 1                     | 0 or 1 (tie-breaking) |
| Weakly informative (dominant ≈ 0.4) | ≈ 0.68          | 1 (0.4 < 0.5)         | 0 (argmax agrees)     |

The v0 soft Hamming gives fractional contributions that reflect uncertainty. The v1-original threshold overcounts (every uncertain position is a difference). The v1-current argmax undercounts (uncertain positions depend on tie-breaking and are often 0 because both profiles tend to peak at the same state - the one favored by the GTR equilibrium frequencies).

For datasets with many uncertain positions (short sequences, high divergence, deep internal nodes), the three approaches produce different initial branch lengths:

- v0: moderate initial guess scaled by uncertainty
- v1-original: inflated initial guess (every uncertain position adds +1)
- v1-current: deflated initial guess (uncertain positions contribute ~0)

## Impact

### Role of initial guess in v1

v1's `initial_guess_mixed()` sets branch lengths directly. The optimization loop then alternates `update_marginal()` (recomputes messages with current branch lengths) and `run_optimize_mixed()` (Newton/grid optimization per edge). The initial guess determines the branch lengths used for the first `update_marginal()` call.

v1 uses Newton's method, which is a local optimizer sensitive to starting point. A poor initial guess can cause:

- Convergence to a different local optimum (branch length surface can have multiple local optima for long branches or weak signal)
- Slower convergence (more iterations to reach the neighborhood of the global optimum)
- Newton failure triggering grid search fallback (which is slower but more reliable)

This contrasts with v0's Brent method, which is a bracketed global optimizer within `[0, MAX_BRANCH_LENGTH]`. v0's soft Hamming distance serves only as a bracket midpoint hint; the final result is insensitive to its exact value.

### Practical severity

For well-resolved trees with strong phylogenetic signal (most positions have sharp profiles), all three approaches give similar initial guesses and the difference is negligible.

For poorly resolved regions (long branches, limited data, divergent sequences), the initial guess quality matters more. These are the same conditions where the branch length surface is least well-behaved for Newton optimization.

The `optimize` command is also used indirectly by `timetree` through the `PartitionOptimizeOps` trait and shared optimization code.

## Data source divergence

Independent of the counting criterion, v1 uses per-edge messages while v0 uses full marginal profiles:

- v1 `edge.msg_to_parent.dis`: subtree likelihood at the child, propagated to the parent through the transition matrix. Carries information from the child's subtree only.
- v1 `edge.msg_to_child.dis`: outgroup likelihood at the parent, propagated to the child through the transition matrix. Carries information from everything except this child's subtree.
- v0 `node.marginal_outgroup_LH`: outgroup likelihood at the child node (from forward pass). Carries rest-of-tree information.
- v0 `node.marginal_subtree_LH`: subtree likelihood at the child node (from backward pass). Carries child-subtree information.

The v0 profiles are evaluated AT the child node (both represent state distributions at the child). The v1 per-edge messages represent state distributions at opposite ends of the edge (parent end and child end). The soft Hamming dot product has different semantics:

- v0 `dot(outgroup_at_child, subtree_at_child)`: measures agreement between two views of the same node's state. When both agree, the state is well-determined. The dot product approximates the posterior certainty.
- v1 `dot(subtree_at_parent, outgroup_at_child)` (hypothetical): would measure agreement between states at different nodes, reflecting both state uncertainty and actual evolutionary change.

Implementing v0's dot product in v1 would require evaluating both messages at the same node (the child), not at opposite ends of the edge.

## Proposed solutions

### S1: Implement v0's profile dot product for dense initial guess

Add a method to the dense partition that computes the soft Hamming distance matching v0:

1. For each edge, obtain the two profile components at the child node (subtree likelihood and outgroup likelihood, both evaluated at the child)
2. Compute `1 - sum(dot(pp[i], pc[i])) / L_eff` across non-gap positions

This requires access to the full node profiles (not per-edge messages). The profiles are available after `update_marginal()` completes - specifically, the node's `profile.dis` contains the combined posterior. The two components can be recovered from the edge messages and the combined profile.

Alternatively, add a dedicated `dense_hamming_distance()` method that works directly on the node posteriors rather than going through `edge_subs()`.

This is the most faithful port but requires profile decomposition that v1's data structures may not directly expose.

### S2: Compute soft Hamming from per-edge messages

Use v1's existing per-edge messages but apply v0's dot product formula:

```
soft_hamming = 1 - sum(dot(msg_to_parent[i], msg_to_child[i])) / L_eff
```

This preserves v1's data source (per-edge messages) but uses v0's counting criterion (soft dot product instead of hard argmax). The resulting distance has different semantics than v0 (dot product of messages at different ends of the edge) but captures uncertainty in the same way.

This is a pragmatic middle ground: easy to implement, captures the key property (uncertain positions contribute fractionally), but not identical to v0.

### S3: Validate current argmax approach empirically

Compare v0, v1-current (argmax), and S2 (soft Hamming) on representative datasets:

- Measure: iteration count to convergence, final log-likelihood, total branch length
- Datasets: flu/h3n2/20 (small, strong signal), ebola (moderate), dengue/2000 (large, weaker signal)
- If argmax produces equivalent final results with acceptable iteration counts, document as an intentional v1 simplification in `docs/port-intentional-changes/`

### S4: Bypass initial guess entirely

Set all initial branch lengths to a fixed heuristic (e.g., `1 / sequence_length` or the total tree length divided by edge count). This avoids the profile-dependent initial guess entirely. The optimization loop should converge to the same result regardless of starting point, assuming Newton/grid search handles it. This is the simplest approach but relies on Newton/grid search handling arbitrary starting points.

## Testing

### T1: Unit test for dense soft Hamming (validates S1 or S2)

Construct a dense partition with known profiles:

- Case A: sharp profiles, parent and child agree at all positions. Expected: soft Hamming ≈ 0, hard Hamming = 0.
- Case B: sharp profiles, parent and child disagree at 5 of 20 positions. Expected: soft Hamming ≈ 0.25, hard Hamming = 0.25. Both agree.
- Case C: uniform profiles at all positions (4-state alphabet). Expected: soft Hamming ≈ 0.75, hard Hamming = variable (depends on argmax tie-breaking). This is the discriminating case.
- Case D: weakly informative profiles (dominant state ≈ 0.4, others ≈ 0.2). Expected: soft Hamming ≈ 0.68, hard Hamming = 0 (argmax agrees). Another discriminating case.

Assert that the soft Hamming implementation produces the expected fractional values for cases C and D.

### T2: Golden master comparison with v0

Run v0's `optimize_tree_marginal` and v1's optimize loop on the same input tree and alignment. Compare:

- Initial branch lengths (after initial guess, before optimization loop)
- Final branch lengths (after convergence)
- Total log-likelihood
- Iteration count to convergence

Capture v0's per-edge `hamming_distance` values from `optimal_t_compressed` as the golden reference. Compare against v1's per-edge initial guess values.

Use a Python capture script to extract v0's intermediate values:

```python
# Capture v0's per-edge soft Hamming distances
for node in tree.find_clades():
    if node.up is None:
        continue
    pp, pc = tt.marginal_branch_profile(node)
    mult = tt.data.multiplicity(mask=node.mask)
    soft_hamming = 1 - np.sum(mult * np.sum(pp * pc, axis=1)) / np.sum(mult)
    # Record (node.name, soft_hamming)
```

Compare against v1's `edge_subs().len() / effective_length` per edge. Positions where the two diverge indicate uncertain sites where the counting criterion matters.

### T3: Convergence sensitivity test

Run the optimize command with two different initial guess strategies (argmax vs soft Hamming) on the same dataset. Compare:

- Whether both converge to the same final log-likelihood (within tolerance)
- Number of iterations required
- Maximum per-edge branch length difference after convergence

If both converge to the same result, the initial guess choice is a performance concern only. If they converge to different results, the initial guess affects correctness.

### T4: Property test for soft Hamming bounds

For any pair of normalized probability profiles `pp` and `pc` over `n` states:

- `0 <= dot(pp, pc) <= 1` (Cauchy-Schwarz, equality when identical)
- `dot(pp, pc) >= 1/n` when both are valid probability distributions (minimum is two uniform distributions)
- `soft_hamming = 1 - dot(pp, pc)` is in `[0, 1 - 1/n]`
- When one profile is a point mass (sharp), `soft_hamming = 1 - pc[argmax(pp)]`, recovering the v1-original threshold test for the specific case where the threshold is the child's probability at the parent's MAP state

## Affected commands

- `optimize` (directly, via `initial_guess_mixed`)
- `timetree` (uses `PartitionOptimizeOps` trait through shared optimization code)

## Related issues

- [Branch length optimization oscillates without damping](M-optimize-oscillation-no-damping.md) - the initial guess quality interacts with oscillation; a poor initial guess combined with no damping amplifies convergence problems
- [Branch mutations have no unified API across partition types](M-core-branch-mutations-no-unified-api.md) - the `PartitionOptimizeOps` trait introduced in PR #470 that changed the dense counting criterion
- Stale JC69 edge messages (fixed) - `update_marginal` was not re-run after replacing the dummy JC69 with the real GTR, so `initial_guess_mixed` read edge messages computed with the wrong model. Fixed by adding `update_marginal` after GTR replacement, matching the ancestral command pattern.
