# The runs folder grows without retention or quota

No app deletes finished runs by itself, and neither has a limit on the number of runs or on their size. The user deletes runs one at a time or as a selection in the run history, which frees the folder of a desktop run and the memory of a web2 run. The desktop app keeps every run in its runs folder on disk (`<app folder>/runs`), so the folder grows with every run that the user keeps. web2 keeps every run of a tab in the in-memory file system of its back-end worker, so the memory of the tab grows with every run until the user deletes runs or the tab reloads or closes; the history shows the size of each run. One ancestral run of `data/sc2/2844` writes about 200 MB of outputs.

> [!WARNING]
> **Needs review.** The `sc2/2844` figures here come from an earlier version of that dataset: a tree with branch lengths in years and an alignment with 340 sequences without a tip. The dataset now holds 2504 genomes on a tree with branch lengths in substitutions per site (`data/sc2/2844/README.md`). Measure the figures again before relying on them.

## Evidence

- `RunStore` in `packages/app-backend/src/runs/store.rs` creates one folder per run under the runs folder (`/treetime/runs` in web2), and `InputStore` keeps the uploads once by content under `.inputs/<sha256>/` of the runs folder. `RunStore::remove_unstarted()` [packages/app-backend/src/runs/store.rs#L172](../../packages/app-backend/src/runs/store.rs#L172) deletes runs that were created for an upload and never started; `RunManager::open()` calls it once at startup for runs older than one day [packages/app-backend/src/runs/manager.rs#L64](../../packages/app-backend/src/runs/manager.rs#L64), which in web2 is once per tab, and the form deletes the upload run of a draft that the user discards
- `DELETE /api/runs/{id}` (`RunManager::delete()`, [packages/app-backend/src/runs/manager.rs#L254](../../packages/app-backend/src/runs/manager.rs#L254)) removes one ended run and the stored inputs that no other run reads. Only the user calls it: no age, size, or count removes a run
- `RunRecord.size` gives the bytes of a run (its outputs and the inputs that only it reads), which the history shows
- The upload limit bounds the inputs of one run, not the outputs or the number of runs
- Each web2 module instance grows its memory up to the 4 GiB of `wasm32` ([kb/decisions/app-web2-wasm-worker-back-end.md](../decisions/app-web2-wasm-worker-back-end.md)), and the limit is per instance: the back-end worker has one instance, and `JobHost` (`packages/app-web2/src/backend/jobs.ts`) starts one job worker with its own instance per run. No quota bounds the memory of all runs of a tab together; the browser decides how much the tab gets, and when the back-end worker runs out, it aborts and every run of the tab is lost. The capacity of the browsers is measured in [N-web2-supported-dataset-size-undocumented.md](N-web2-supported-dataset-size-undocumented.md)

## Possible changes

- A quota on the total size of the runs: new runs are refused with a clear error when it is full
- Retention that deletes finished runs older than a set age, except pinned runs

> [!IMPORTANT]
> **Decision required.** No bound beyond the user's own deletes is chosen. Open choices: a quota, an age-based retention, or both; whether either applies to desktop, whose users keep their runs, or only to the memory of a web2 tab; and whether the UI warns before a run is removed. Evidence: one large run writes about 200 MB, and no other bound limits the runs.
