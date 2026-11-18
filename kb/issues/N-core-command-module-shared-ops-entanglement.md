# Command orchestration mixes application policy, domain workflow, and I/O

> [!IMPORTANT]
> **Decision required.** The application operation boundary is undecided: either the shared application crate `app-commands` owns validated operations and returns in-memory results, or each adapter invokes domain pipelines directly through transport-neutral request types. See "Open design question".

The command boundary is not a thin adapter. Command runners combine argument translation, input loading, pipeline invocation, output policy, serialization, and filesystem writes, while one timetree pipeline function sequences the complete scientific workflow.

## Evidence

- `fn run_ancestral_reconstruction()` reads FASTA, maps CLI arguments, runs inference, builds output projections, and writes several formats [packages/app-commands/src/commands/ancestral/run.rs#L48-L185](../../packages/app-commands/src/commands/ancestral/run.rs#L48-L185).
- The same application-level shape appears in clock, homoplasy, mugration, optimize, prune, and timetree runners under [`packages/app-commands/src/commands`](../../packages/app-commands/src/commands).
- `fn timetree::pipeline::run()` [packages/treetime/src/timetree/pipeline.rs#L59-L146](../../packages/treetime/src/timetree/pipeline.rs#L59-L146) sequences the complete scientific workflow through step functions: date loading, clock estimation, the pre-loop steps (ML branch-length optimization, rerooting, clock filter), coalescent initialization, refinement, confidence intervals, and result assembly.

- `app-commands` is already the shared application crate: `app-cli`, `app-server`, `app-backend`, and `app-napi` depend on it, and `app-wasm` reaches it through `app-backend`. The adapters run a command through `CommandArgs::execute()` [packages/app-commands/src/command.rs#L135-L145](../../packages/app-commands/src/command.rs#L135-L145), which returns `()`: the runners write their results to files, and the adapters do not receive them in memory. The one exception is `clock` in the CLI, which calls `run_clock()` directly to print the regression chart from its result [packages/app-cli/src/run.rs#L114-L120](../../packages/app-cli/src/run.rs#L114-L120).

Domain modules do not import `commands/`, but the boundary between application orchestration and file output inside `app-commands` is not defined.

## Open design question

The application operation boundary must serve CLI, HTTP, N-API, desktop, and future Python clients. Moving all orchestration into the CLI crate would leave the other clients without an owner. The unresolved choice is whether `app-commands` owns validated operations and returns in-memory results, or each adapter invokes domain pipelines directly through transport-neutral request types.

Moving `commands/` waits until this boundary is decided. Scientific stage ordering and fallback behavior must remain unchanged unless separately approved.

## Related issues

- [N-command-output-ownership-is-scattered.md](N-command-output-ownership-is-scattered.md)

## Proposed implementation plan

All cited locations are current. The plan stops at the open decision; implement one option after approval. Both options keep the stage order and fallbacks of `timetree::pipeline::run()` unchanged, and both start with the same step.

1. Common step: split each runner in `packages/app-commands/src/commands/*/run.rs` into an operation that takes parsed inputs and returns a result struct, and a writer that turns the result into files. `run_clock()` already returns a result and serves as the pattern; coordinate the writer half with [N-command-output-ownership-is-scattered.md](N-command-output-ownership-is-scattered.md), which owns output policy
2. Option A, `app-commands` owns operations: make `CommandArgs::execute()` return a result enum with one variant per command, move file writing behind an output request that the adapter passes in, and let `app-cli` print the clock chart from the returned result instead of calling `run_clock()` directly. `app-server`, `app-backend`, and `app-napi` then return results in memory without reading files back
3. Option B, adapters call domain pipelines: define transport-neutral request types in `app-commands` (or a new contracts crate), let each adapter call the domain pipeline (`timetree::pipeline::run()` and the siblings) and its own output step, and reduce `app-commands` to request parsing and validation
4. Tests: per command, one test of the operation half on a small in-memory dataset that asserts the result struct (oracle: the files that the current runner writes for the same input, compared after the split as a golden master), and one wiring test per adapter entry point
5. Validate: `just test-rs`, `just lint-rs`, `just test-wasm`, `just test-wasm-mt`, and `just smoke`, because the split must not change any output file
