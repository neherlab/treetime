# Dissolve hacks/ module

## Description

`hacks/` contains one 14-line function `fix_branch_length` that clamps branch lengths to a minimum. Used by 3 callers: `partition/marginal_passes.rs`, `partition/marginal_dense.rs`, `gtr/infer_gtr/common.rs`. Module name normalizes technical debt.

## Fix

Move `fix_branch_length` to `seq/alignment.rs` or inline at the 3 call sites. Delete the `hacks/` module entirely.

## Related issues

- Source: [kb/issues/M-core-remaining-architectural-debt-after-extraction.md](../issues/M-core-remaining-architectural-debt-after-extraction.md)
- Also tracked in: [C3 zero-branch-length-clamping](convention-zero-branch-length-clamping.md)
