# Reset form keeps the uploads of the pending run

The apps upload the input files of a new analysis into a run that is created and not started (status `created`), the pending upload run of the draft. "Reset form" replaces the draft with a fresh one and drops `upload_run_id`, but releases none of the uploads of that run. The next upload creates a new pending run, and the old run keeps its files.

In the server and the desktop app, the files stay on disk until the run manager next opens the runs folder and removes unstarted runs older than one day. In web2, the runs folder lives in the memory of the back-end worker for the life of the tab, so every reset after a large upload keeps that upload in memory until the page reloads.

## Evidence

- `reset()` of the draft store replaces the draft and clears the selections without releasing uploads [packages/app-ui/src/store/draft.ts#L84-L91](../../packages/app-ui/src/store/draft.ts#L84-L91)
- The upload code releases an upload only while the draft still names its run: `releaseFromPendingRun()` returns when `upload_run_id` differs [packages/app-ui/src/analysis/useRunUploads.ts#L42-L49](../../packages/app-ui/src/analysis/useRunUploads.ts#L42-L49). An example config that is still loading when the user resets the form leaves its uploaded files in the old run for the same reason
- `RunManager::open()` removes unstarted runs created more than `UNSTARTED_RUN_LIFETIME` (one day) before it opens [packages/app-backend/src/runs/manager.rs#L41](../../packages/app-backend/src/runs/manager.rs#L41), [#L56-L58](../../packages/app-backend/src/runs/manager.rs#L56-L58)
- web2 keeps runs in the memory of the tab: [kb/decisions/app-web2-wasm-worker-back-end.md](../decisions/app-web2-wasm-worker-back-end.md)

> [!IMPORTANT]
> **Investigation required.** Measure the memory of the web2 back-end worker after several resets that each follow an upload of a large alignment (for example `sc2/5100`), to find whether the kept uploads can approach the 4 GiB limit of the module.

## Fix direction

Release the uploads of the pending run when the form resets (or delete the unstarted run), in the draft reset path, so every app frees them at once.

## Validation

- A UI test resets a draft that holds uploads and checks that the back end receives a release of each upload of the old run
