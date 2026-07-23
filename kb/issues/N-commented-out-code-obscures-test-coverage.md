# Commented-out code obscures test coverage and supported behavior

Maintained Rust sources contain commented-out parameterized test cases and a commented-out `Seq::splice()` implementation [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L17`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L17) [`packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L18`](../../packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L18) [`packages/treetime-primitives/src/seq.rs#L197`](../../packages/treetime-primitives/src/seq.rs#L197).

## Failure mode

- Commented `#[case]` rows look like coverage but are never collected, reported, or counted as ignored by the test runner.
- A reader cannot tell whether each case is unsupported, slow, flaky, redundant, or temporarily failing.
- Commented production functions can silently become incompatible with the surrounding types while appearing to document intended capability.

## Required disposition

Every disabled case needs an explicit outcome: restore it as an executable case, mark it ignored with a concrete reason when the runner must retain it, link it to a known issue when it exposes a defect, or delete it. Delete the commented `Seq::splice()` implementation unless an active issue specifies its required contract.

## Validation

- Test listing confirms restored and ignored cases are collected under descriptive names.
- Each retained ignore has a self-contained reason and issue link where applicable.
- Repository search finds no commented test attributes or commented function bodies in maintained Rust sources.
