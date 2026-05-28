# Ancestral parallel traversals repeat identical error-capture boilerplate

Three parallel traversal sites in ancestral reconstruction repeat a 10-line error-capture pattern: create `Arc<Mutex<Option<Report>>>`, execute `par_iter` closure with error capture, call `extract_parallel_error`.

## Locations

- `packages/treetime/src/ancestral/marginal.rs:83:` forward pass error capture
- `packages/treetime/src/ancestral/marginal.rs:119:` backward pass error capture
- `packages/treetime/src/ancestral/fitch.rs:201:` Fitch pass error capture

All three follow the identical pattern:

```rust
let error: Arc<Mutex<Option<Report>>> = Arc::new(Mutex::new(None));
graph.par_iter_...(|mut node| {
  if let Err(e) = operation(&mut node) {
    let mut guard = error.lock();
    if guard.is_none() { *guard = Some(e); }
    return GraphTraversalContinuation::Stop;
  }
  GraphTraversalContinuation::Continue
});
extract_parallel_error(error)?;
```

## Impact

Pure maintainability. The pattern is correct and consistent. Adding a new parallel traversal requires copying the boilerplate correctly.

## Action

Extract `par_traversal_with_error(graph, direction, op) -> Result<()>` helper that encapsulates the `Arc<Mutex<Option<Report>>>` + closure + `extract_parallel_error` pattern.

## Related

- `kb/issues/M-timetree-parallel-traversal-unwrap-panic.md` -- timetree traversals use `.unwrap()` instead of this pattern (policy divergence)
