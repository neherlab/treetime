# Initial branch length guess always overwrites input values

`initial_guess_mixed()` at [packages/treetime/src/commands/optimize/optimize_unified.rs#L300-L328](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L300-L328) runs unconditionally before the Newton iteration loop. It computes `edge_subs().len() / edge_effective_length()` for each edge and overwrites the branch length, discarding whatever values the input tree provided.

## Problem

The `#subs/length` formula is a crude estimate based on discrete mutation counts. For trees with well-calibrated input branch lengths (e.g., from RAxML, IQ-TREE, or a previous TreeTime run), the input values are a better starting point for Newton iteration than the discrete-count estimate. Overwriting them forces Newton to converge from a worse initial position.

The formula is also sensitive to ambiguous and gap-heavy positions. Although gap positions are excluded via `edge_effective_length()`, the discrete mutation count from `edge_subs()` can differ from the continuous ML estimate, producing initial values that are systematically biased relative to the input.

## Expected behavior

The initial guess should be optional. When the input tree has branch lengths, the optimizer should use them as the starting point for Newton iteration. The initial guess should run only when branch lengths are missing or explicitly requested.

## Current call site

`run_optimize()` at [packages/treetime/src/commands/optimize/run.rs#L131](../../packages/treetime/src/commands/optimize/run.rs#L131):

```rust
initial_guess_mixed(&graph, &mixed_partitions)?;
```

Called unconditionally after marginal reconstruction and before the optimization loop.

## Proposed solution

Add a CLI flag `--initial-guess` (default: `auto`) with three modes:

- `auto`: run `initial_guess_mixed()` only when any edge has `branch_length == None`
- `always`: run unconditionally (current behavior)
- `never`: skip, use input branch lengths as-is

The `auto` mode checks whether all edges have `Some(branch_length)`. If so, the input values are used directly. If any edge lacks a branch length, the initial guess fills in the missing values.

## Impact

Low. The optimizer converges to the same ML branch lengths regardless of starting point, given enough iterations. The difference is convergence speed: starting from input branch lengths requires fewer Newton steps than starting from the discrete-count estimate.

## v0 handling

v0's `optimize_tree()` does not have a separate initial guess step. It starts Newton/Brent optimization from the current branch lengths, which are either input values or values from a previous optimization round.
