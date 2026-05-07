# Remove or implement build_branch_distributions todo stub

`build_branch_distributions()` in [packages/treetime/src/commands/timetree/inference/branch_distributions.rs](../../packages/treetime/src/commands/timetree/inference/branch_distributions.rs) is exported as a public API but its body is `todo!()`. The function includes a `Joint` variant from the `BranchDistributionMethod` enum that was removed in v1 (joint reconstruction was intentionally dropped).

## Impact

Any caller of `build_branch_distributions()` panics at runtime. The function is not called from the current timetree pipeline (which uses `create_poisson_branch_distributions` instead), so the panic is not reachable through normal operation. The dead public surface creates confusion about the API contract.

## Affected code

- Stub: [packages/treetime/src/commands/timetree/inference/branch_distributions.rs](../../packages/treetime/src/commands/timetree/inference/branch_distributions.rs)
- The `Joint` variant of `BranchDistributionMethod` enum in the same file references the removed joint reconstruction mode

## Fix

Either remove `build_branch_distributions` and the `Joint` variant entirely, or replace the `todo!()` with the actual implementation if the function is needed. If kept, remove the `Joint` variant to match the v1 design decision documented in [decisions/ancestral-joint-reconstruction-removed.md](../decisions/ancestral-joint-reconstruction-removed.md).

## Related issues

- Source: [M-timetree-branch-distributions-todo-stub.md](../issues/M-timetree-branch-distributions-todo-stub.md) -- delete after full resolution
- [ancestral-joint-reconstruction-removed.md](../decisions/ancestral-joint-reconstruction-removed.md) -- joint reconstruction removal
- [N-timetree-n-branches-posterior-unimplemented.md](../issues/N-timetree-n-branches-posterior-unimplemented.md) -- related unimplemented branch distribution feature
