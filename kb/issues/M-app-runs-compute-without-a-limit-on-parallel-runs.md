# Runs compute without a limit on parallel runs

Both apps start every run at once, and nothing limits how many runs compute at the same time. In web2 each run computes in a job worker of its own, with the threaded module at 8 threads by default and its own linear memory of up to 4 GiB, so a few large runs started together multiply the threads and the memory that the tab asks of the browser. The desktop app computes every run in a thread of its own, and all runs share one pool of processing threads, so parallel runs compete for the cores and add up their memory.

## Evidence

- web2: `JobHost` in `packages/app-web2/src/backend/jobs.ts` starts one job worker for every start order of the back end, with no queue and no limit; each job worker creates its own instance of a module ([kb/decisions/app-web2-wasm-worker-back-end.md](../decisions/app-web2-wasm-worker-back-end.md), "One job worker per run"). The thread count of each run is the URL parameter `jobs`, default 8 (`DEFAULT_JOBS` in `packages/app-web2/src/backend/threads.ts`)
- Desktop: `ThreadExecutor` in `packages/app-backend/src/bridge/executor.rs` computes each run in a thread of its own; the threads share the global rayon pool
- `RunList.active_runs` in `packages/app-backend/src/runs/record.rs` reports the number of runs computing now, but no code compares it with a limit

## Possible changes

- A configurable number of runs computing at once, with further started runs waiting in a queue
- The UI shows a waiting run as queued, with its position

> [!IMPORTANT]
> **Decision required.** Open choices, each independent:
>
> - The default number of runs computing at once per app: no limit (today), a fixed number, or a number from the core count (desktop) and the thread count of a run (web2). Evidence: a web2 run takes 8 threads and its own memory of up to 4 GiB; desktop runs share one rayon pool
> - A waiting run: cancellable like a running one, or only removable by delete
> - Threads in web2: each run keeps the `jobs` value of the page, or the page divides its threads between the runs that compute at the same time
