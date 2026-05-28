# Extract parallel traversal error-capture helper from ancestral boilerplate

Three ancestral traversal sites duplicate the 10-line `Arc<Mutex<Option<Report>>>` error-capture pattern.

## Current state

`ancestral/marginal.rs:83,119` and `ancestral/fitch.rs:201` each manually create the error arc, wrap the closure, and call `extract_parallel_error`.

## Target state

A `par_traversal_with_error(graph, direction, op) -> Result<()>` helper encapsulates the pattern. All three ancestral sites and the two timetree traversal sites (currently using `.unwrap()`) use the helper.

## Implementation

1. Define `par_traversal_with_error` in a shared location (e.g., `ancestral/` or a graph utilities module)
2. The helper takes a graph reference, traversal direction (forward/backward), and a fallible operation closure
3. Internally creates the `Arc<Mutex<Option<Report>>>`, executes the parallel traversal, and calls `extract_parallel_error`
4. Replace the 3 ancestral sites with calls to the helper
5. Optionally: convert the 2 timetree traversal sites (`forward_pass.rs:17`, `backward_pass.rs:26`) from `.unwrap()` to the helper (see related ticket `timetree-traversals-unwrap-on-fallible-distribution-math.md`)

## Related issues

Source: `kb/issues/N-ancestral-parallel-traversal-error-boilerplate.md`
Related: `kb/issues/M-timetree-parallel-traversal-unwrap-panic.md`
