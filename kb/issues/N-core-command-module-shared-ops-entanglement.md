# Command orchestration mixes application policy, domain workflow, and I/O

> [!IMPORTANT]
> **Decision required.** The application operation boundary is undecided: either the shared application crate `app-commands` owns validated operations and returns in-memory results, or each adapter invokes domain pipelines directly through transport-neutral request types. See "Open design question".

The command boundary is not a thin adapter. Command runners combine argument translation, input loading, pipeline invocation, output policy, serialization, and filesystem writes, while one timetree pipeline function sequences the complete scientific workflow.

## Evidence

- `fn run_ancestral_reconstruction()` reads FASTA, maps CLI arguments, runs inference, builds output projections, and writes several formats [packages/app-commands/src/commands/ancestral/run.rs#L42-L160](../../packages/app-commands/src/commands/ancestral/run.rs#L42-L160).
- The same application-level shape appears in clock, mugration, optimize, prune, and timetree runners under [`packages/app-commands/src/commands`](../../packages/app-commands/src/commands).
- `fn timetree::pipeline::run()` [packages/treetime/src/timetree/pipeline.rs#L54-L124](../../packages/treetime/src/timetree/pipeline.rs#L54-L124) sequences the complete scientific workflow through step functions: date loading, clock estimation, the pre-loop steps (ML branch-length optimization, rerooting, clock filter), coalescent initialization, refinement, confidence intervals, and result assembly.

- `app-commands` is already the shared application crate: `app-cli`, `app-server`, `app-backend`, and `app-napi` depend on it, and `app-wasm` reaches it through `app-backend`. The adapters run a command through `CommandArgs::execute()` [packages/app-commands/src/command.rs#L133-L143](../../packages/app-commands/src/command.rs#L133-L143), which returns `()`: the runners write their results to files, and the adapters do not receive them in memory. The one exception is `clock` in the CLI, which calls `run_clock()` directly to print the regression chart from its result [packages/app-cli/src/run.rs#L48-L54](../../packages/app-cli/src/run.rs#L48-L54).

Domain modules do not import `commands/`, but the boundary between application orchestration and file output inside `app-commands` is not defined.

## Open design question

The application operation boundary must serve CLI, HTTP, N-API, desktop, and future Python clients. Moving all orchestration into the CLI crate would leave the other clients without an owner. The unresolved choice is whether `app-commands` owns validated operations and returns in-memory results, or each adapter invokes domain pipelines directly through transport-neutral request types.

Moving `commands/` waits until this boundary is decided. Scientific stage ordering and fallback behavior must remain unchanged unless separately approved.

## Related issues

- [N-command-output-ownership-is-scattered.md](N-command-output-ownership-is-scattered.md)
