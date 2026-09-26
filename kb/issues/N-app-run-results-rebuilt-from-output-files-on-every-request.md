# Run results and the Auspice document are rebuilt from the output files on every request

`GET /api/runs/{id}/results` and `GET /api/runs/{id}/auspice` read the output files of a finished run and build their answer again on every request, although the files of a finished run do not change. On large datasets each request costs more than a tenth of a second of back-end time, and the comparison of two runs pays it twice.

## Evidence

- `fn run_results()` in `packages/app-commands/src/results/run_results.rs` reads the run record, then `fn results_of_record()` parses the Auspice JSON, the tree and the other output files of `out/` and assembles `RunResults`
- `fn run_auspice()` in `packages/app-commands/src/results/auspice.rs` parses the Auspice JSON with `fn read_auspice()` (`packages/app-commands/src/results/outputs.rs`) and computes the color scales with `fn display_auspice()`
- `fn compare_runs()` in `packages/app-commands/src/results/compare.rs` calls `run_results` for both runs, and `fn clade_in_runs()` in `packages/app-commands/src/results/clades.rs` reads the result tree of every other finished time-tree run on each request
- Nothing caches these answers in the back end; `AppService` (`packages/app-commands/src/bridge/service.rs`) passes each call straight to these functions

## Measurements

Release N-API addon, requests sent through an in-process router, median of repeated calls. The time is spent in Rust; the transport adds less than 1 %.

| Request                             | Answer size | Time in Rust |
| ----------------------------------- | ----------- | ------------ |
| `run_results`, `flu/h3n2/500`       | 336 KB      | 31.5 ms      |
| `run_auspice`, `flu/h3n2/500`       | 167 KB      | 26.1 ms      |
| `run_results`, `dengue/2000`        | 1.65 MB     | 151 ms       |
| `run_auspice`, `dengue/2000`        | 1.10 MB     | 129 ms       |

## Impact

- Every page that opens the results of a large run waits this long for each answer, and every invalidation of the run's cache entries (for example a title change) reads the files again
- Comparing two large runs and searching a clade across many runs scale with the number and size of the runs read

## Possible changes

- Keep the built answers of finished runs in a bounded in-memory cache keyed by run id, dropped when the run is deleted or purged
- Write the answers next to the outputs when a run ends, and serve the stored JSON

The choice between the two is open. A cache keeps the answers consistent with changes to the result-building code; stored answers survive a restart but go stale when that code changes.
