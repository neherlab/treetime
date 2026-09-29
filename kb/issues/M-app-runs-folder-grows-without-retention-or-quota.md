# The runs folder grows without retention or quota

The server never deletes runs and has no limit on the number of runs or on the size of the runs folder. One ancestral run of `data/sc2/2844` writes about 200 MB of outputs. On the nightly server, anyone with the password can fill the disk by starting runs in a loop, and the disk also fills over time without abuse.

## Evidence

- `RunStore` in `packages/app-commands/src/runs/store.rs` creates one folder per run under `--runs-dir` and has no removal by age or size
- The upload limit (`--max-upload-size`, 50 MB per run on the nightly server) bounds the inputs of one run, not the outputs or the number of runs
- `docs/dev/developer_guide.md`, section "Nightly web deploy", states that the server never deletes runs

## Possible changes

- A quota on the total size of the runs folder: new runs are refused with a clear error when it is full
- Retention that deletes finished runs older than a set age, except pinned runs

Both must be off by default, because local and desktop users keep their runs, and on only in the deploy. Open design questions: default age and quota for the deploy, and whether the UI warns before a run is removed.
