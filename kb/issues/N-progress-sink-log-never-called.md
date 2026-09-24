# `ProgressSink::log` is never called

## Summary

`ProgressSink::log` in `packages/treetime/src/progress.rs` is a public trait method that no workspace code calls. Its only caller is the `progress_log!` macro in the same file, and no crate uses that macro. The CLI, NAPI, and server sinks implement `log`, but library pipelines never send log messages through the sink.

## Consequence

`just check-all` fails on a clean tree: the `pub_unused_in_workspace` dylint lint reports the method as unused, and the gate denies warnings.

## Open question

Choose one of these directions:

- **Route pipeline logging through the sink**: call `progress_log!` from the library pipelines, so the desktop and web clients receive log messages the same way they receive progress
- **Remove the unused logging path**: delete `log`, `log_enabled`, `progress_log!`, and the sink implementations of these methods
