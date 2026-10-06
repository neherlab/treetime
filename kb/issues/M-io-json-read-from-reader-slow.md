# JSON files are read through a slow streaming parser

`json_read()` in `packages/treetime-utils/src/io/json.rs` parses JSON with `serde_json::Deserializer::from_reader()` wrapped in `serde_stacker`. The streaming reader parses large files 8 to 20 times slower than reading the file into memory and parsing the byte slice. The app results endpoint (`GET /api/runs/{id}/results`) reads output files through `json_read_file()`, so large runs take tens of seconds to show their results.

## Evidence

Dev profile, test build, files written by app runs with default settings:

| File | `json_read_file()` | `fs::read()` + `serde_json::from_slice()` |
| --- | --- | --- |
| `homoplasy.stats.json` of `data/sc2/4500` (123 MB) | 18.1 s | 0.93 s |
| `homoplasy.stats.json` of `data/mpox/clade-ii/2000` (221 MB) | 32.9 s | 4.1 s |
| `homoplasy.auspice.json` of `data/sc2/4500` | 6.2 s | not measured |
| `homoplasy.auspice.json` of `data/mpox/clade-ii/2000` | 3.3 s | not measured |

Both reads deserialize into the same read type of the homoplasy results (`struct HomoplasyStatsFile` in `packages/app-commands/src/results/homoplasy.rs`), which skips the ranked list of ambiguous changes. Building the results from the parsed file takes less than 0.1 s. `GET /api/runs/{id}/results` takes 20.5 s for the `sc2/4500` homoplasy run and 30.8 s for the `mpox/clade-ii/2000` run; every results view reads its Auspice file the same way, so ancestral, timetree, and the other views of large runs are slow for the same reason.

## Fix direction

Read the decompressed file into memory inside `json_read()` and parse it with `serde_json::from_slice()`, keeping the `serde_stacker` wrapper for deeply nested trees. The memory cost is one copy of the decompressed file during parsing.

> [!IMPORTANT]
> **Investigation required.** Measure the change on the release profile and on deeply nested Auspice trees, and check the peak memory of the largest outputs, before changing the shared reader.
