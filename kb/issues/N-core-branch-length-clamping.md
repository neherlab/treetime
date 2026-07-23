# Marginal branch-length clamping has unresolved numerical policy and ownership

`fn fix_branch_length()` raises every branch below `MIN_BRANCH_LENGTH_FRACTION / sequence_length` to that threshold before marginal message propagation [`packages/treetime/src/hacks/fix_branch_length.rs#L4`](../../packages/treetime/src/hacks/fix_branch_length.rs#L4). The name does not reveal that it clamps, the public `hacks` namespace does not identify the owning algorithm, and the helper is used only by marginal passes [`packages/treetime/src/partition/marginal_passes.rs#L148`](../../packages/treetime/src/partition/marginal_passes.rs#L148) [`packages/treetime/src/partition/marginal_passes.rs#L398`](../../packages/treetime/src/partition/marginal_passes.rs#L398).

The clamp avoids non-finite values in the current log-space path at zero branch length. Mathematically, $e^{Q0}=I$ is valid, so changing the clamp is a numerical and scientific decision rather than architecture cleanup.

## Separable work

- Preserve behavior while moving the helper beside marginal inference as a private, policy-revealing operation.
- Decide clamping, exact-zero handling, or rejection only with parity and numerical validation.

## Related tickets

- [kb/tickets/architecture-dissolve-hacks-module.md](../tickets/architecture-dissolve-hacks-module.md)
- [kb/tickets/convention-zero-branch-length-clamping.md](../tickets/convention-zero-branch-length-clamping.md)
