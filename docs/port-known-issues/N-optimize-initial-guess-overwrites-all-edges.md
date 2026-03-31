# Auto mode initial guess overwrites all edges

When `--initial-guess=auto` detects any edge with a missing branch length, `initial_guess_mixed()` overwrites ALL edges including those with valid input branch lengths. A tree with 99 well-calibrated edges and 1 missing edge has all 100 overwritten with the discrete-count estimate.

## Expected behavior

`initial_guess_mixed()` should only fill in edges with missing branch lengths, preserving valid input values on other edges.

## Current behavior

`initial_guess_mixed()` at [packages/treetime/src/commands/optimize/optimize_unified.rs#L300-L328](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L300-L328) iterates over all edges unconditionally.

## Impact

Negligible. The optimizer converges to the same ML branch lengths regardless of starting point. Mixed trees (some valid, some missing) are rare in practice. Starting from the discrete-count estimate adds a few extra Newton iterations compared to starting from the input values.

## Proposed fix

Add a predicate inside `initial_guess_mixed()` to skip edges that already have a valid (non-None, non-NaN) branch length.
