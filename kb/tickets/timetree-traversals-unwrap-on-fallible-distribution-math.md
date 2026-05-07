# Backward and forward traversals use unwrap on fallible distribution math

The backward and forward marginal passes in the timetree inference pipeline use `.unwrap()` inside parallel traversal closures. Distribution operations (convolution, multiplication, division) return `Result`, but the traversal closures discard the error via `.unwrap()`, converting recoverable errors into panics.

## Affected code

- Backward pass: `.unwrap()` on distribution math results inside `backward_pass` traversal closure
- Forward pass: `.unwrap()` on distribution math results inside `forward_pass` traversal closure
- `compute_node_contributions`: `.unwrap()` inside breadth-first traversal closure at [packages/treetime/src/commands/timetree/coalescent/contributions.rs#L62-L70](../../packages/treetime/src/commands/timetree/coalescent/contributions.rs#L62-L70)

## Impact

When a distribution operation fails (degenerate input, numerical overflow, empty distribution), the entire process panics instead of producing a diagnostic error. Since these closures run inside parallel traversals (`par_iter_breadth_first_forward/backward`), the panic aborts the rayon thread pool.

Triggering conditions include:

- Trees with very short or zero branch lengths (producing degenerate `expQt` matrices)
- Numerical underflow after many narrow distribution multiplications
- `Tc` evaluation failure in coalescent contribution computation

## Fix

Replace `.unwrap()` calls with error collection: either propagate errors out of the traversal by collecting into a `Vec<Report>` and returning them after the traversal completes, or refactor the traversal API to support fallible callbacks. The `compute_node_contributions` case at [contributions.rs#L62-L70](../../packages/treetime/src/commands/timetree/coalescent/contributions.rs#L62-L70) uses `.unwrap()` because `iter_breadth_first_forward` expects an infallible callback.

## Related issues

- Source: [M-timetree-inference-unwrap-in-traversals.md](../issues/M-timetree-inference-unwrap-in-traversals.md) -- delete after full resolution
- [M-optimize-negative-branch-length-validation.md](M-optimize-negative-branch-length-validation.md): missing negative branch-length validation, a trigger for degenerate distribution math
