# Runs compute without a limit on parallel runs

The server starts every run at once, and all runs share one pool of processing threads (`-j`). Nothing limits how many runs compute at the same time. A few large runs started together can exceed the memory of the machine, and the kernel then stops the whole app, including every other run. This applies to `treetime-server` on a self-hosted machine and to the desktop app. The hosted web app runs the static web2 build, in which each run computes in its own WebAssembly worker.

## Evidence

- `AppService::start_run()` in `packages/app-backend/src/bridge/service.rs` hands each run to `RunExecutor::execute()`; `ThreadExecutor` spawns one thread per run with `thread::Builder` [packages/app-backend/src/bridge/executor.rs#L22-L29](../../packages/app-backend/src/bridge/executor.rs#L22-L29); the threads share the global rayon pool set by `--jobs`
- `RunList.active_runs` in `packages/app-backend/src/runs/record.rs` (filled by `active_count()` at [packages/app-backend/src/runs/manager.rs#L88](../../packages/app-backend/src/runs/manager.rs#L88)) reports the number of runs computing now, but no code compares it with a limit
- Request rate limits do not bound this cost: one request starts a run of any size, and a limit strict enough to protect the server blocks normal use

## Possible changes

- A configurable number of runs computing at once, with further started runs waiting in a queue, and a maximum number of waiting runs answered with HTTP 429 when full
- The UI shows a waiting run as queued, with its position

Open design questions: queue size, whether a waiting run can be cancelled like a running one, and whether the desktop app uses the same limit.

## Proposed implementation plan

No code limits the number of computing runs. The open design questions above (queue size, cancel of a waiting run, desktop limit) stay open; the plan works for every answer:

1. Put the limit into the owner of run execution: a queueing `RunExecutor` in `packages/app-backend/src/bridge/` that wraps `ThreadExecutor`, holds started runs in a FIFO queue, and starts the next one when a run thread ends (the run thread signals its end on every outcome, panic included). Configure the limit and the queue size through the server settings
2. Reject a start with HTTP 429 when the queue is full, through the existing error mapping of `packages/app-server`
3. Add a `queued` run status with its position to the run record and the contracts (`just gen`), and show it in the run views of `packages/app-ui`
4. Tests in `packages/app-backend/src/bridge/__tests__/`: with a limit of 1 and gated runs (as `GatedReader` gates uploads), the second run stays queued until the first ends, and a full queue rejects the next start. Oracle: the configured limit. See red: a limit that the executor ignores
5. Validate: `just lint-rs`, `just test-rs`, `just test-ts`, `just generated-check`
