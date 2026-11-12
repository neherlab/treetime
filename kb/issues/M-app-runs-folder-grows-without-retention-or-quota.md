# The runs folder grows without retention or quota

No app deletes finished runs, and neither has a limit on the number of runs or on their size. The desktop app keeps every run in its runs folder on disk (`<app folder>/runs`), so the folder grows with every run until the user deletes run folders by hand. web2 keeps every run of a tab in the in-memory file system of its back-end worker, so the memory of the tab grows with every run until the tab reloads or closes, and no run can be freed without losing the others. One ancestral run of `data/sc2/2844` writes about 200 MB of outputs.

> [!WARNING]
> **Needs review.** The `sc2/2844` figures here come from an earlier version of that dataset: a tree with branch lengths in years and an alignment with 340 sequences without a tip. The dataset now holds 2504 genomes on a tree with branch lengths in substitutions per site (`data/sc2/2844/README.md`). Measure the figures again before relying on them.

## Evidence

- `RunStore` in `packages/app-backend/src/runs/store.rs` creates one folder per run under the runs folder (`/treetime/runs` in web2). Its only removal is `RunStore::remove_unstarted()` [packages/app-backend/src/runs/store.rs#L133-L144](../../packages/app-backend/src/runs/store.rs#L133-L144), which deletes runs that were created for an upload and never started. `RunManager::open()` calls it once at startup for runs older than one day [packages/app-backend/src/runs/manager.rs#L67-L69](../../packages/app-backend/src/runs/manager.rs#L67-L69), which in web2 is once per tab. Runs that started are never removed, by age or by size
- No route, `RunManager`, or `RunStore` method deletes a run
- The upload limit bounds the inputs of one run, not the outputs or the number of runs
- A web2 module instance grows its memory up to the 4 GiB of `wasm32` ([kb/decisions/app-web2-wasm-worker-back-end.md](../decisions/app-web2-wasm-worker-back-end.md)), so the runs of one tab together stay below that bound or the back-end worker aborts and every run of the tab is lost

## Possible changes

- A delete route for single runs, so a user frees the disk space of a desktop run and the memory of a web2 run without losing the other runs
- A quota on the total size of the runs: new runs are refused with a clear error when it is full
- Retention that deletes finished runs older than a set age, except pinned runs

> [!IMPORTANT]
> **Decision required.** Beyond a delete route, no bound is chosen. Open choices: a quota, an age-based retention, or both; whether either applies to desktop, whose users keep their runs, or only to the memory of a web2 tab; and whether the UI warns before a run is removed. Evidence: one large run writes about 200 MB, and no other bound limits the runs.
