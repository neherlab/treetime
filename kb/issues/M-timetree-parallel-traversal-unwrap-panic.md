# Timetree parallel traversals use unwrap instead of error capture

Timetree forward and backward distribution propagation use `.unwrap()` inside parallel traversal closures, while every other parallel traversal in the codebase (5 sites in ancestral) uses the error-capture pattern (`Arc<Mutex<Option<Report>>>` + `extract_parallel_error`).

A panic inside a rayon parallel closure aborts the process.

## Locations

- `packages/treetime/src/timetree/inference/forward_pass.rs:17:` `propagate_distributions_forward_single_node(&mut node).unwrap()`
- `packages/treetime/src/timetree/inference/backward_pass.rs:26:` `propagate_distributions_backward_single_node(&mut node, coalescent_contributions).unwrap()`

## Expected behavior

Use the same error-capture pattern as ancestral traversals:

```rust
let error: Arc<Mutex<Option<Report>>> = Arc::new(Mutex::new(None));
graph.par_iter_breadth_first_forward(|mut node| {
  if let Err(e) = propagate_distributions_forward_single_node(&mut node) {
    let mut guard = error.lock();
    if guard.is_none() { *guard = Some(e); }
    return GraphTraversalContinuation::Stop;
  }
  GraphTraversalContinuation::Continue
});
extract_parallel_error(error)?;
```

## Reference

Ancestral error-capture sites: `packages/treetime/src/ancestral/marginal.rs:98,120,142,178:` and `packages/treetime/src/ancestral/fitch.rs:196:`

## Related

Duplicate detection report `.memory/08-dupe-scan/synthesis1.md` finding #3 notes this as a policy divergence (A10).
