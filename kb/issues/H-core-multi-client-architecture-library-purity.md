# Application command boundary does not isolate clients

TreeTime has CLI, HTTP, N-API, desktop, and web clients, but no application boundary owns command requests, operation execution, cancellation, and in-memory results for all of them.

## Evidence

- `app-api` forwards six command functions directly to `treetime::commands` and re-exports concrete pipeline and progress types [`packages/app-api/src/commands.rs#L14`](../../packages/app-api/src/commands.rs#L14) [`packages/app-api/src/pipelines.rs#L1`](../../packages/app-api/src/pipelines.rs#L1).
- `app-cli` calls `treetime::commands` directly, while `app-server` and `app-napi` call through `app-api`; the purported boundary is therefore optional rather than authoritative [`packages/app-cli/src/bin/treetime.rs#L23`](../../packages/app-cli/src/bin/treetime.rs#L23) [`packages/app-server/src/routes.rs#L42`](../../packages/app-server/src/routes.rs#L42) [`packages/app-napi/src/commands.rs#L27`](../../packages/app-napi/src/commands.rs#L27).
- `app-server` still imports concrete model, reroot, optimization, and command types from `treetime` to construct command arguments [`packages/app-server/src/args.rs#L4`](../../packages/app-server/src/args.rs#L4).
- The command runners are filesystem-oriented and return concrete core types, so callers inherit CLI-era input and output decisions instead of consuming a stable application request/result contract.

## Impact

- Adding or changing one command requires coordinated edits in core commands, `app-api`, server routes, N-API tasks, transport schemas, and client dispatch.
- Cancellation, progress, validation, and error semantics vary by adapter.
- Moving command orchestration into the CLI crate would deepen the mismatch because non-CLI clients also need the same operations.

## Design axes

### Application ownership

- O1. Make a shared application crate own transport-neutral requests, operation execution, progress/cancellation, and in-memory results. CLI, HTTP, and N-API become adapters.
- O2. Remove the pass-through `app-api` surface and let every adapter call domain pipelines directly. Shared transport contracts remain separate, and cross-adapter parity must be enforced by tests.

### Core `clap` dependency

Keeping `clap` derives in core types is an existing project decision. This issue does not reopen that choice; it concerns operation ownership and client isolation.

No implementation ticket is ready until application ownership is selected.

## Related issues

- [H-core-command-module-shared-ops-entanglement.md](H-core-command-module-shared-ops-entanglement.md)
- [H-app-transport-contracts-diverge-across-clients.md](H-app-transport-contracts-diverge-across-clients.md)
- [H-app-napi-cancellation-is-process-global.md](H-app-napi-cancellation-is-process-global.md)
