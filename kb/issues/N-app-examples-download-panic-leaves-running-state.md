# A panic in the example download leaves the download state at running

`ExampleDownloads::start()` in `packages/app-backend/src/examples_download.rs` starts a thread that downloads and unpacks the example datasets and publishes the result as the download state, which the app polls. The thread's `JoinHandle` is dropped, so nothing observes the thread when it ends.

`download_examples()` returns every expected failure (network, file system, archive) as an error, and the thread publishes it. A panic inside the thread ends it before it publishes anything, so the state stays at `running`, and the app shows a download that never ends. The app must restart before it can start another download.

The run threads of `packages/app-backend/src/bridge/service.rs` do not have this gap: `run_job()` in `packages/app-backend/src/job.rs` catches a panic of the run with `catch_unwind` and records it as the outcome of the run.

## Fix direction

Catch a panic of the download work in the thread, as `run_job()` does, and publish it as a failed download. The `spawn_handle_dropped` expectation on `ExampleDownloads::start()` then keeps its reason without the gap.

## Proposed implementation plan

`fn ExampleDownloads.start()` ([packages/app-backend/src/examples_download.rs#L75-L98](../../packages/app-backend/src/examples_download.rs#L75-L98)) still spawns `fn ExampleDownloads.run()` with the handle dropped, and `run()` ([packages/app-backend/src/examples_download.rs#L100-L124](../../packages/app-backend/src/examples_download.rs#L100-L124)) publishes only the `Ok`/`Err` of `download_examples()`. The panic catch of the run threads now sits in `fn run_job_step()` ([packages/app-backend/src/job.rs#L36-L46](../../packages/app-backend/src/job.rs#L36-L46)).

1. In `fn ExampleDownloads.run()`, wrap the `download_examples(...)` call in `catch_unwind(AssertUnwindSafe(...))`, and map a panic payload to `ExamplesDownload::Failed` with the text of `treetime_utils::error::panic_message()` ([packages/treetime-utils/src/error.rs#L4](../../packages/treetime-utils/src/error.rs#L4)), so every end of the thread publishes a terminal state
2. Make the download function injectable for the test only if `ExampleDownloads` cannot take a replacement download function yet (for example through a constructor that takes it); keep the production path unchanged otherwise
3. Test in `packages/app-backend/src/__tests__/`: a download function that panics; wait on the published `ExamplesDownload::Failed` event of `runs.app_events()` (no sleep), and assert the whole `ExamplesDownloadStatus`. Oracle: the contract above (every ended download leaves `Running`). See red: remove the `catch_unwind` and confirm that the state stays `Running`
4. Keep the `spawn_handle_dropped` expectation and extend its reason with the panic catch

### Plan validation

- `just lint-rs` and `just test-rs -E 'package(app-backend)'` pass
