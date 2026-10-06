# JSON files are read through a slow streaming parser

`json_read_file()` in `packages/treetime-utils/src/io/json.rs` parses JSON with `serde_json::Deserializer::from_reader()` over the `Box<dyn BufRead>` that `read_file_with()` passes to it. This path executes about 150 CPU instructions per input byte in the release profile and about 760 in the dev profile. Reading the file into memory and parsing it with `serde_json::from_slice()` executes 12 to 18 instructions per byte in the release profile and about 50 in the dev profile. The app results endpoint (`GET /api/runs/{id}/results`) reads the output files of a run through `json_read_file()`, so on the dev server a large run takes 10 to 30 seconds to show its results.

## Cause

Three causes multiply:

- **The stream reader of serde_json reads one byte per call.** `IoRead` pulls each byte through `std::io::Bytes`, updates a line and column counter for each byte, and copies each string byte into a scratch buffer. The slice reader scans string content 8 bytes at a time, borrows strings that have no escapes, and computes the line and column only for an error message. The documentation of `serde_json::from_reader()` states that reading the whole file into memory and parsing it with `from_slice()` is usually faster
- **The boxed reader disables the single-byte fast path of std.** `open_file_or_stdin()` returns `Box<dyn BufRead>` because the input can be a file or standard input. std reads one byte from a `BufReader<R>` with an inlined buffer check, but this specialization does not apply through `Box<dyn BufRead>`. Each byte is then a virtual `read()` call with a 1-byte buffer, through `BufReader::read()`, `Buffer::fill_buf()`, and `<&[u8] as Read>::read()`. Iterating the bytes of the 125 MB statistics file of `sc2/4500` through this reader, without parsing, takes 9.6 s in the dev profile
- **The dev profile does not optimize this path.** `[profile.dev.package."*"] opt-level = 2` applies to the code that a dependency compiles itself. The generic and `#[inline]` functions of std and serde_json are compiled in the workspace crate that calls them, at opt-level 0 and with the debug precondition checks of `slice::from_raw_parts()`. In a perf profile of the dev build, the 1-byte read calls take about 85% of the samples and the JSON parser about 10%

The writers do not have the second cause: `create_file_or_stdout()` keeps the concrete `BufWriter<Compressor>` outermost, so the boxed writer inside it is called once per 256 KiB buffer flush.

The `serde_stacker` wrapper adds no measurable time, and it is necessary. Without it, the dev build overflows the 2 MiB stack of a `spawn_blocking` thread on the Auspice tree of `sc2/4500`, which is 157 nodes deep.

## Evidence

Standalone benchmark on a thread with a 2 MiB stack, with the files in the page cache. Each file is deserialized into the read type of the results page (`HomoplasyStatsFile` in `packages/app-commands/src/results/homoplasy.rs`, `AuspiceTree` in `packages/treetime-io/src/auspice_types.rs`). The columns are the current `json_read_file()`, the same parser over a concrete `BufReader<Decompressor>`, and the decompressed file read into memory and parsed with `from_slice()` behind `serde_stacker`.

Instructions per input byte, release profile:

| File                                                    | Current | `BufReader<Decompressor>` | In memory |
| ------------------------------------------------------- | ------- | ------------------------- | --------- |
| `homoplasy.stats.json` of `sc2/4500` (125 MB)           | 154     | 82                        | 12        |
| `homoplasy.stats.json` of `mpox/clade-ii/2000` (227 MB) | 154     | 81                        | 14        |
| Auspice JSON of `sc2/4500` (48 MB)                      | 149     | 76                        | 13        |
| `timetree.auspice.json` of `mpox/clade-ii/500` (199 MB) | 154     | 81                        | 18        |

Time, minimum of repeated runs on a loaded host:

| File                            | Dev: current | Dev: `BufReader` | Dev: in memory | Release: current | Release: `BufReader` | Release: in memory |
| ------------------------------- | ------------ | ---------------- | -------------- | ---------------- | -------------------- | ------------------ |
| `sc2/4500` statistics           | 9.8 s        | 4.6 s            | 0.63 s         | 1.17 s           | 0.68 s               | 0.20 s             |
| `mpox/clade-ii/2000` statistics | 17.7 s       | 8.6 s            | 1.45 s         | 3.45 s           | 1.86 s               | 0.48 s             |
| `sc2/4500` Auspice              | 3.7 s        | 1.6 s            | 0.41 s         | 0.73 s           | 0.37 s               | 0.12 s             |
| `mpox/clade-ii/500` Auspice     | 15.7 s       | 7.5 s            | 1.56 s         | 2.48 s           | 1.34 s               | 0.40 s             |

Peak resident memory, release profile:

| File                            | Current | In memory, buffer grown by `read_to_end()` | In memory, buffer sized from the file length |
| ------------------------------- | ------- | ------------------------------------------ | -------------------------------------------- |
| `mpox/clade-ii/2000` statistics | 11 MB   | 272 MB                                     | 232 MB                                       |
| `mpox/clade-ii/500` Auspice     | 100 MB  | 361 MB                                     | 293 MB                                       |
| `sc2/4500` Auspice              | 30 MB   | 95 MB                                      | 76 MB                                        |

The statistics files are large because of the ranked list of ambiguous changes, which the results page skips ([N-homoplasy-stats-json-size.md](N-homoplasy-stats-json-size.md)).

## Fix direction

> [!IMPORTANT]
> **Decision required.** Two independent changes are possible, alone or together:
>
> - **Parse from memory**: read the decompressed input in `json_read_file()` into a buffer and parse it with `serde_json::from_slice()` behind `serde_stacker`. This is 6 to 16 times faster than the current path. The peak memory grows by the size of the decompressed file while the parser runs, for example from 11 MB to 232 MB for the `mpox/clade-ii/2000` statistics. A buffer sized from the file length avoids the larger peak of a growing buffer, but compressed files and standard input have no known length
> - **Concrete reader**: return a concrete `BufReader<Decompressor>` from `open_file_or_stdin()` for both files and standard input, with standard input read through a `Decompressor` without compression. The boxed reader inside `Decompressor` is then called once per 256 KiB buffer refill instead of once per byte. This halves the time of the streaming parser without a memory cost and applies to every reader of `read_file_with()`, but the streaming parser stays 3 to 7 times slower than parsing from memory
