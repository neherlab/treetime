# Indel rate estimation returns zero on all-zero-BL input trees

On an all-zero branch length input tree with identical sequences and indels present, the optimization loop can converge before the Poisson indel term is ever active. The indel-bearing edge ends at a substitution-driven value near `min_branch_length` instead of the Poisson MLE $t^* = k / \mu$.

## Trigger conditions

All four must hold simultaneously:

1. All input branch lengths are exactly zero
2. Sequences are identical (no substitution signal)
3. Indels are present on at least one edge
4. `InitialGuessMode::Auto` (default) preserves the zero BLs

## Location

- [packages/treetime/src/commands/optimize/optimize_indel.rs#L20-L23](../../packages/treetime/src/commands/optimize/optimize_indel.rs#L20-L23): `poisson_indel_log_lh()` early return when `mu == 0.0` - silently drops the indel term
- [packages/treetime/src/commands/optimize/optimize_unified.rs#L289](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L289): `estimate_indel_rate()` called once per `run_optimize_mixed` invocation - returns 0 when total BL is 0
- [packages/treetime/src/commands/optimize/run.rs#L132-L133](../../packages/treetime/src/commands/optimize/run.rs#L132-L133): `Auto` mode passes `overwrite_valid=false`, preserving zero BLs
- [packages/treetime/src/commands/optimize/run.rs#L163](../../packages/treetime/src/commands/optimize/run.rs#L163): convergence check uses substitution-only likelihood from `update_marginal`

## Mechanism

### Iteration 1

1. `update_marginal` at BL=0 → substitution log-likelihood $L_1$
2. `|L_1 - \text{MIN}| > dp$ → continue
3. `run_optimize_mixed`:
   - `estimate_indel_rate()` → total BL is 0 → returns 0.0
   - `poisson_indel_log_lh(k, 0.0, t)` returns neutral metrics → indel term inactive
   - Newton optimizes substitution-only. For identical sequences, optimum is at zero. Newton pushes toward zero, clamped at `min_branch_length = one_mutation * 0.01`
4. `apply_damping(0.75, i=0)`: BL ≈ `one_mutation * 0.01 * 0.25`

### Iteration 2

1. `update_marginal` at BL ≈ $7.8 \times 10^{-5}$ → substitution log-likelihood $L_2$
2. For JC69 with identical sequences: $\Delta \log L \approx 16 \times 7.8 \times 10^{-5} \approx 0.00125 < dp = 0.01$
3. Convergence check fires. Loop exits. `run_optimize_mixed` never runs with positive `indel_rate`.

### Result

The indel-bearing edge has BL ≈ $7.8 \times 10^{-5}$ (substitution-driven) instead of $t^* = k / \mu$ (indel-driven). The Poisson term was never active during optimization.

## Impact

Negligible. The four trigger conditions are individually uncommon and almost never co-occur in practice. Input trees from phylogenetic inference tools (RAxML, IQ-TREE, FastTree) have positive branch lengths, making conditions 1 and 2 inapplicable. The result is a tiny positive BL rather than zero, so the tree topology is preserved.

## Fix options

- F1. Treat exact zero BLs as invalid in `InitialGuessMode::Auto` when any edge carries indels, so `initial_guess_mixed()` bootstraps positive BLs before `run_optimize_mixed()`
- F2. Make `poisson_indel_log_lh(k > 0, mu = 0, t)` signal an invalid state (return $-\infty$ or error) instead of neutral metrics, forcing a rate bootstrap
- F3. Include the indel contribution in the outer-loop convergence check alongside substitution likelihood
- F4. Run a rate-bootstrap pass before the main loop: compute initial BLs, estimate indel rate, then enter the convergence loop with a positive rate

## v0 comparison

v0 does not model indels in branch length optimization. Not applicable.
