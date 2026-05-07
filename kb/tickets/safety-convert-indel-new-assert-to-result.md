# Convert InDel::new assert! to Result

## Description

`packages/treetime/src/seq/indel.rs:23-24:`

`assert!(range.0 <= range.1)` and `assert_eq!(seq.len(), range.1 - range.0)` panic on invalid input instead of returning `Result`. Reachable from any command that processes alignments with indels when input data contains malformed gap annotations.

Any malformed input reaching `InDel::new` causes an unrecoverable panic. Production users see panic backtraces instead of actionable error messages.

## Fix

Replace `assert!` calls in `InDel::new` with `Result` return type and proper error reporting.

## Related issues

- Source: [N-production-unwrap-expect-audit.md](../issues/N-production-unwrap-expect-audit.md) -- delete after full resolution
