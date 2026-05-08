# Move From<BranchSplitArgs> impl from clock/ to commands/clock/

## Description

`clock/find_best_root/params.rs` imports `commands::clock::args::BranchSplitArgs` for a `From<&BranchSplitArgs> for BranchPointOptimizationParams` impl. Domain module depends on CLI args struct.

## Fix

Move the `From` impl to `commands/clock/args.rs` or a new `commands/clock/params_from_args.rs`. The domain module `clock/` should define `BranchPointOptimizationParams` without knowing about CLI types. The command layer converts.

One file change, one import removal.

## Related issues

- Source: [M-core-remaining-architectural-debt-after-extraction](../issues/M-core-remaining-architectural-debt-after-extraction.md)
- Also tracked in: [M-core-cli-library-separation-blockers](../issues/M-core-cli-library-separation-blockers.md) (B1 subset)
