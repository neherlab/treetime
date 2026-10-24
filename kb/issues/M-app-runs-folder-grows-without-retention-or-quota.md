# The runs folder grows without retention or quota

The server never deletes finished runs and has no limit on the number of runs or on the size of the runs folder. One ancestral run of `data/sc2/2844` writes about 200 MB of outputs. On the nightly server, anyone with the password can fill the disk by starting runs in a loop, and the disk also fills over time without abuse.

> [!WARNING]
> **Needs review.** The `sc2/2844` figures here come from an earlier version of that dataset: a tree with branch lengths in years and an alignment with 340 sequences without a tip. The dataset now holds 2504 genomes on a tree with branch lengths in substitutions per site (`data/sc2/2844/README.md`). Measure the figures again before relying on them.

## Evidence

- `RunStore` in `packages/app-backend/src/runs/store.rs` creates one folder per run under `--runs-dir`. Its only removal is `RunStore::remove_unstarted()` [packages/app-backend/src/runs/store.rs#L133-L144](../../packages/app-backend/src/runs/store.rs#L133-L144), which deletes runs that were created for an upload and never started. `RunManager::open()` calls it once at startup for runs older than one day [packages/app-backend/src/runs/manager.rs#L67-L69](../../packages/app-backend/src/runs/manager.rs#L67-L69). Runs that started are never removed, by age or by size
- The upload limit (`--max-upload-size`, 50 MB per run on the nightly server) bounds the inputs of one run, not the outputs or the number of runs
- `docs/dev/developer_guide.md`, section "Nightly web deploy", states that the server never deletes runs

## Possible changes

- A quota on the total size of the runs folder: new runs are refused with a clear error when it is full
- Retention that deletes finished runs older than a set age, except pinned runs

Both must be off by default, because local and desktop users keep their runs, and on only in the deploy.

> [!IMPORTANT]
> **Decision required.** The limits for the deploy are not chosen. Open choices: a quota, an age-based retention, or both; the default age and quota for the deploy; and whether the UI warns before a run is removed. Evidence: one large run writes about 200 MB, and the server has no other bound on the runs folder.
