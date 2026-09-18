# Comparison harness builds the measured revision with its own tooling

The output-equality harness checks out the `rust` baseline into a worktree and builds it by running that revision's `./dev/dev` and `./dev/docker/run`. A revision under measurement therefore defines the mechanism intended to contain its own build execution.

`fn build_baseline_binary()` adds a detached worktree at the baseline commit and invokes the checkout's own build wrapper [dev/compare-baseline#L288-L296](../../dev/compare-baseline#L288-L296). This makes the measured revision part of the build and containment mechanism.

## Decision axes

### Containment owner

- O1. Use a separately identified trusted harness revision to build the baseline binary.
- O2. Execute the baseline revision's own build tooling. This matches how that revision builds but grants the measured revision control of containment.

Recommendation: O1. The baseline commit's build behavior can be inspected separately without executing it.

## Recommendation

Build the baseline binary from a trusted harness revision rather than the checked-out baseline's own scripts. Keep this issue ticketless until the containment-owner policy is approved.
