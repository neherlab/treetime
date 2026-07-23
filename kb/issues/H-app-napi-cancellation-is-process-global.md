# N-API cancellation is process-global instead of task-local

Every N-API command task shares one static cancellation flag [`packages/app-napi/src/progress.rs#L8`](../../packages/app-napi/src/progress.rs#L8). Each progress sink reads that flag, and every macro-generated task resets it before starting [`packages/app-napi/src/progress.rs#L63`](../../packages/app-napi/src/progress.rs#L63) [`packages/app-napi/src/commands.rs#L64`](../../packages/app-napi/src/commands.rs#L64).

## Data flow

The exported `cancel()` operation writes one `CANCELLED: AtomicBool`. All six command task families construct `NapiProgressSink` values that read the same bit through `is_cancelled()`. Starting any macro-generated task calls `reset_cancel()` against that shared bit. `AncestralTaskNoop` follows a separate compute path and does not perform the same reset.

The atomic protects concurrent access to one value; it does not associate the value with a task identity or prevent one task from changing another task's control state.

## Failure mode

- Cancelling one task cancels every concurrent task.
- Starting a new task clears a cancellation request made for an older task.
- The no-op task path does not follow the same reset sequence.

## Required contract

Cancellation ownership must identify one job. Construct a task-local cancellation token, retain access to that token in the corresponding exported task handle or job registry, and inject it into that task's progress sink. Task startup must never mutate cancellation state owned by another task.

The transport contract must define:

- how jobs are identified;
- whether cancellation is idempotent;
- the terminal result reported for a cancelled job;
- what happens when cancellation arrives before start or after completion.

## Validation

- Run two concurrent tasks, cancel one, and prove the other continues.
- Cancel an older task while a new task starts and prove the request is not cleared.
- Exercise every generated command family and the no-op task through the same lifecycle contract.
- Verify success, failure, panic, and cancellation each produce exactly one terminal result.

## Related issues

- [H-app-transport-contracts-diverge-across-clients.md](H-app-transport-contracts-diverge-across-clients.md)
- [H-core-multi-client-architecture-library-purity.md](H-core-multi-client-architecture-library-purity.md)
- [kb/tickets/app-unify-web-desktop-and-typescript-transport-contracts.md](../tickets/app-unify-web-desktop-and-typescript-transport-contracts.md)
