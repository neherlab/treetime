# Standalone helper optimizers duplicate production logic and drift from unified behavior

`run_optimize_dense()` and `run_optimize_sparse()` duplicate the branch-length optimization logic implemented by the unified optimizer, but they are only used by tests. This duplication already diverged in scientifically important ways, especially in the zero-branch shortcut, and it makes dense-sparse equivalence tests an unreliable proxy for the production `optimize` path.

## Problem

The repository currently has three optimizer implementations for branch lengths:

- unified production path in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L180-L250`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L180-L250)
- standalone dense helper in [`packages/treetime/src/commands/optimize/optimize_dense.rs#L77-L165`](../../packages/treetime/src/commands/optimize/optimize_dense.rs#L77-L165)
- standalone sparse helper in [`packages/treetime/src/commands/optimize/optimize_sparse.rs#L142-L224`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L142-L224)

The helper implementations duplicate:

- zero-branch shortcut logic
- Newton update logic
- fallback grid-search logic
- branch-length update structure

The duplication is not harmless. The zero-branch logic already drifted:

- unified path uses `is_zero_branch_optimal()` in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L153-L173`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L153-L173)
- standalone dense reimplements a similar thresholded shortcut inline in [`packages/treetime/src/commands/optimize/optimize_dense.rs#L104-L123`](../../packages/treetime/src/commands/optimize/optimize_dense.rs#L104-L123)
- standalone sparse reimplements a different shortcut in [`packages/treetime/src/commands/optimize/optimize_sparse.rs#L160-L175`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L160-L175), with a different threshold and no derivative check

This means a future bug fix can easily land in the production optimizer while leaving the helper optimizers behind.

## Blast radius

### Affected code

- [`packages/treetime/src/commands/optimize/optimize_unified.rs#L180-L250`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L180-L250)
- [`packages/treetime/src/commands/optimize/optimize_dense.rs#L77-L165`](../../packages/treetime/src/commands/optimize/optimize_dense.rs#L77-L165)
- [`packages/treetime/src/commands/optimize/optimize_sparse.rs#L142-L224`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs#L142-L224)

### Caller status

Exhaustive search on `dev` found no production callers of `run_optimize_dense()` or `run_optimize_sparse()`. They are used only by dense-sparse equivalence tests:

- [`packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs)
- [`packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_convergence.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_convergence.rs)
- [`packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_validity.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_validity.rs)

### User-visible effect

The duplication does not currently affect CLI behavior directly, because the production `optimize` command uses the unified path. The damage is indirect:

- tests can disagree with production behavior for reasons unrelated to dense-vs-sparse science
- duplicated bugs can survive in helper code and confuse future debugging
- fixes to convergence or zero-branch handling need to be repeated in multiple places

## Related known issues

- [Zero-branch shortcut underflows and depends on a fixed likelihood threshold](M-optimize-zero-branch-shortcut-raw-likelihood-threshold.md)
- [Standalone sparse optimizer skips derivative check for zero branches](N-optimize-sparse-zero-branch-no-derivative.md)

This issue is about the duplication itself. The two linked issues are concrete manifestations of that duplication.

## Proposed solutions

### S1. Share the zero-branch decision logic

Extract one shared helper for zero-branch detection and reuse it from unified, dense, and sparse paths.

Pros:

- removes the most obvious existing divergence
- keeps behavior aligned across production and helper code
- reduces the chance of fixing [the zero-branch shortcut issue](M-optimize-zero-branch-shortcut-raw-likelihood-threshold.md) in one place only

Cons:

- still leaves Newton and grid-search logic duplicated

### S2. Move helper optimizers behind `#[cfg(test)]` and minimize their surface

Keep the helper optimizers only for tests and mark them clearly as non-production utilities.

Pros:

- matches current caller reality
- reduces confusion about which optimizer is authoritative

Cons:

- does not by itself remove duplicated behavior

### S3. Delete the helper optimizers and route tests through the unified optimizer

Rewrite the dense-sparse equivalence tests to use the unified optimizer with single-type partition vectors, then remove the helper implementations.

Pros:

- one source of truth for optimization behavior
- tests measure the actual production path

Cons:

- tests lose one axis of direct helper-vs-helper comparison
- requires test refactoring

### Recommended direction

S1 is the minimum correction. S3 is the cleanest end state. If the helper optimizers stay, they should not carry their own distinct zero-branch rules.

## Proposed tests and verification

### T1. Shared zero-branch behavior

If a shared helper is introduced, add tests that feed the same synthetic coefficients through unified, dense, and sparse wrappers and assert the same zero-branch decision.

### T2. Production-helper consistency

For a small dense-only or sparse-only setup, compare the branch lengths produced by the standalone helper against the unified optimizer configured with the same single partition type.

### T3. Test intent clarity

Update dense-sparse equivalence tests so they state whether they validate:

- scientific equivalence of dense and sparse contributions
- or exact optimizer-path equivalence

That distinction matters once helper code and production code no longer drift silently.
