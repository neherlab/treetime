# Dissolve hacks/ module

## Description

`hacks/` contains `fix_branch_length`, which applies the marginal-inference branch-length floor. The current callers are the dense and sparse marginal passes.

## Fix

Move the helper beside the marginal passes, give it a name that states the floor policy, and delete the `hacks/` module. Preserve the existing numerical behavior; changing the floor or its call sites requires a separate scientific decision.

## Related issues

- Source: [kb/issues/N-core-branch-length-clamping.md](../issues/N-core-branch-length-clamping.md)
- Related: [kb/tickets/convention-zero-branch-length-clamping.md](convention-zero-branch-length-clamping.md)
