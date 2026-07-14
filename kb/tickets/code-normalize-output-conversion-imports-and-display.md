# Normalize output conversion imports and display

Apply the repository’s Rust conventions to the output conversion surface as one cleanup task.

## Required changes

- Move command conversion imports from six `From::from` bodies in `packages/app-server/src/args.rs` to module scope.
- Replace handwritten `Display` implementations for `OutputSelection` and `TopologyOrderTargetSourceArg` with `strum::Display` derives and explicit kebab-case formatting.
- Preserve the exact Clap and diagnostic spellings for every variant.

## Validation

- Parameterized display/Clap spelling tests for every variant.
- Server argument conversion tests.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/N-code-quality-conventions.md](../issues/N-code-quality-conventions.md)
