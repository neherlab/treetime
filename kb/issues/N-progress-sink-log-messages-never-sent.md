# Progress sinks never receive log messages

`trait ProgressSink` declares `log()` and `log_enabled()` next to `report()` [`packages/treetime/src/progress.rs#L5-L9`](../../packages/treetime/src/progress.rs#L5-L9), and the `progress_log!`, `progress_error!`, `progress_warn!`, `progress_info!`, `progress_debug!` and `progress_trace!` macros route a message to `log()` when `log_enabled()` allows it [`packages/treetime/src/progress.rs#L38-L80`](../../packages/treetime/src/progress.rs#L38-L80). No code invokes any of these macros, and nothing else calls `log()`.

Four sinks implement `log()` for their front ends:

- `BarProgress` and `TextProgress` for the command-line interface [`packages/app-cli/src/cli/progress.rs#L36-L108`](../../packages/app-cli/src/cli/progress.rs#L36-L108)
- `ChannelProgress` for the HTTP server's event stream [`packages/app-server/src/sse.rs#L133-L152`](../../packages/app-server/src/sse.rs#L133-L152)
- `NapiProgressSink` for the Node addon [`packages/app-napi/src/progress.rs#L39-L58`](../../packages/app-napi/src/progress.rs#L39-L58)

Those implementations never run, so a front end that relies on the progress sink for diagnostics receives only stage and fraction updates.

## Decision axes

- O1. Send core diagnostics through the sink: replace the relevant `log` crate calls in core operations with the `progress_*!` macros, so every front end receives the same messages
- O2. Remove `log()`, `log_enabled()`, the macros and the four implementations, and keep diagnostics on the `log` crate only
