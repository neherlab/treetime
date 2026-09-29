# Runs compute without a limit on parallel runs

The server starts every run at once, and all runs share one pool of processing threads (`-j`). Nothing limits how many runs compute at the same time. On the nightly server (2 CPUs, 3 GB memory limit in `dev/deploy/hetzner/compose.yaml`), a few large runs started together can exceed the memory limit, and the kernel then stops the whole app, including every other run.

## Evidence

- `AppService::start_run()` in `packages/app-commands/src/bridge/service.rs` spawns one thread per run with `thread::Builder`; the threads share the global rayon pool set by `--jobs`
- `RunList.active_runs` in `packages/app-commands/src/runs/record.rs` reports the number of runs computing now, but no code compares it with a limit
- Request rate limits do not bound this cost: one request starts a run of any size, and a limit strict enough to protect the server blocks normal use

## Possible changes

- A configurable number of runs computing at once (1 on the nightly server), with further started runs waiting in a queue, and a maximum number of waiting runs answered with HTTP 429 when full
- The UI shows a waiting run as queued, with its position

Open design questions: queue size, whether a waiting run can be cancelled like a running one, and whether the desktop app uses the same limit.
