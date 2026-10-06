# JSON files are read through a slow streaming parser

`json_read_file()` in `packages/treetime-utils/src/io/json.rs` parses JSON with `serde_json::Deserializer::from_reader()` over the `Box<dyn BufRead>` that `read_file_with()` passes to it. This reader does a large amount of work for each input byte, and the output files that the apps read are large because they are indented ([M-io-json-output-indentation-dominates-file-size.md](M-io-json-output-indentation-dominates-file-size.md)). On the dev server, the results page of a `homoplasy` run of `sc2/4500` waits 15.0 s for `GET /api/runs/{id}/results`, then 4.3 s for `GET /api/runs/{id}/auspice`. A perf profile of the server during the first request attributes at least 70% of the samples to parsing the statistics file with the stream reader.

## Cause

- **The stream reader of serde_json handles one byte per call.** `IoRead` gets each byte through `peek()` or `next()`, which call `LineColIterator::next()` (a line and column counter) and `std::io::Bytes::next()`, which returns an `Option<io::Result<u8>>`. Whitespace, string content, and skipped values go through this chain one byte at a time. The slice reader of serde_json skips whitespace with direct indexing, scans string content 8 bytes at a time, and borrows strings that have no escapes. serde_json's `Read` trait is sealed, so a caller cannot supply a reader that scans its buffer in bulk
- **The boxed reader adds a virtual call per byte.** std reads one byte from a `BufReader<R>` with an inlined buffer check (rust-lang/rust#116785), but only when the type is exactly `BufReader<R>`. Through `Box<dyn BufRead>`, each byte is a virtual `read()` call with a 1-byte buffer, through `BufReader::read()`, `Buffer::fill_buf()`, and `<&[u8] as Read>::read()`
- **The build profile multiplies the cost.** In the dev profile, the generic and `#[inline]` functions of std and serde_json are compiled in the workspace crate that calls them, at opt-level 0 and with the debug precondition checks of `slice::from_raw_parts()`; `[profile.dev.package."*"] opt-level = 2` does not apply to them. In the release profile (incremental, 16 codegen units), the per-byte functions stay separate calls. In the `dist` profile (fat LTO), they are inlined, and whitespace skipping (`parse_whitespace()`) takes 58% of the time for an indented Auspice file

These are not causes:

- **`Decompressor` on an uncompressed file.** `CompressionType::None` wraps the `BufReader<File>` in a `Box<dyn Read>`. The outer 256 KiB `BufReader` calls it once per refill, and the inner `BufReader` passes reads of its own size or larger straight to the file. A 47.8 MB file takes 183 `read()` system calls of 256 KiB each
- **`serde_stacker`.** It adds no measurable time, and it is necessary: without it, the dev build overflows the 2 MiB stack of a `spawn_blocking` thread on the Auspice tree of `sc2/4500`, which is 157 nodes deep

## Evidence

Benchmark binary on a thread with a 2 MiB stack, files in the page cache, deserialized into the read types of the results page (`HomoplasyStatsFile` in `packages/app-commands/src/results/homoplasy.rs`, `AuspiceTree` in `packages/treetime-io/src/auspice_types.rs`). Three readers:

- **Current**: `json_read_file()`
- **Concrete**: the same parser over a concrete `BufReader<Decompressor>`
- **In memory**: the decompressed file read into a buffer and parsed with `from_slice()` behind `serde_stacker`

Files: `homoplasy.stats.json` of `sc2/4500` and of `mpox/clade-ii/2000`, the Auspice JSON of `sc2/4500`, and a `timetree.auspice.json` of `mpox/clade-ii/500`. "Indented" is the file as written, "compact" is the same JSON without whitespace, and "compact, no ranked list" also drops `ambiguous.all.ranked`, which the results page skips ([N-homoplasy-stats-json-size.md](N-homoplasy-stats-json-size.md)).

Time in the dev profile, which the dev server uses:

| File                                                     | Size    | Current | Concrete | In memory |
| -------------------------------------------------------- | ------- | ------- | -------- | --------- |
| `sc2/4500` statistics, indented                          | 125 MB  | 9.8 s   | 4.6 s    | 0.63 s    |
| `sc2/4500` statistics, compact                           | 77 MB   | 6.2 s   | 2.4 s    | 0.23 s    |
| `sc2/4500` statistics, compact, no ranked list           | 4.5 MB  | 0.42 s  | 0.20 s   | 0.06 s    |
| `mpox/clade-ii/2000` statistics, indented                | 227 MB  | 17.7 s  | 8.6 s    | 1.45 s    |
| `mpox/clade-ii/2000` statistics, compact, no ranked list | 7.6 MB  | 0.74 s  | 0.34 s   | 0.12 s    |
| `sc2/4500` Auspice, indented                             | 47.8 MB | 3.7 s   | 1.6 s    | 0.41 s    |
| `sc2/4500` Auspice, compact                              | 1.3 MB  | 0.26 s  | 0.10 s   | 0.06 s    |
| `mpox/clade-ii/500` Auspice, indented                    | 199 MB  | 15.7 s  | 7.5 s    | 1.56 s    |
| `mpox/clade-ii/500` Auspice, compact                     | 17.8 MB | 1.83 s  | 0.88 s   | 0.29 s    |

Time in the `profiling` profile (the `dist` settings with debug information):

| File                                           | Current | Concrete | In memory |
| ---------------------------------------------- | ------- | -------- | --------- |
| `sc2/4500` statistics, indented                | 0.96 s  | 0.34 s   | 0.13 s    |
| `sc2/4500` statistics, compact, no ranked list | 0.05 s  | 0.03 s   | 0.02 s    |
| `mpox/clade-ii/2000` statistics, indented      | 1.82 s  | 1.03 s   | 0.43 s    |
| `sc2/4500` Auspice, indented                   | 0.56 s  | 0.17 s   | 0.11 s    |
| `mpox/clade-ii/500` Auspice, indented          | 1.38 s  | 0.48 s   | 0.30 s    |
| `mpox/clade-ii/500` Auspice, compact           | 0.24 s  | 0.14 s   | 0.11 s    |

The times are the minimum of repeated runs on a shared host. CPU instructions per input byte are stable across files: about 150 for the current reader, 80 for the concrete reader, and 12 to 18 for the in-memory reader in the release profile; 129, 44, and 17 in the `profiling` profile.

Peak resident memory, release profile:

| File                                      | Current | In memory, buffer grown by `read_to_end()` | In memory, buffer sized from the file length |
| ----------------------------------------- | ------- | ------------------------------------------ | -------------------------------------------- |
| `mpox/clade-ii/2000` statistics, indented | 11 MB   | 272 MB                                     | 232 MB                                       |
| `mpox/clade-ii/500` Auspice, indented     | 100 MB  | 361 MB                                     | 293 MB                                       |
| `sc2/4500` Auspice, indented              | 30 MB   | 95 MB                                      | 76 MB                                        |

A serde_json build with the opt-in `BufferedIoRead` of serde-rs/json#1294, which runs the slice parser over a 256 KiB buffer, streams as fast as `from_slice()` at opt-level 3 (`sc2/4500` statistics, indented: 0.11 s against 0.13 s; `mpox/clade-ii/500` Auspice, indented: 0.25 s against 0.27 s) and about 2 times slower than `from_slice()` at opt-level 0. Removing the line and column counter alone makes the concrete reader 17 to 37% faster at opt-level 3, and does not change it at opt-level 0.

## Fix direction

The largest change is outside the reader: compact output files and a statistics file without the ranked list of ambiguous changes reduce the bytes to parse by 10 to 37 times ([M-io-json-output-indentation-dominates-file-size.md](M-io-json-output-indentation-dominates-file-size.md), [N-homoplasy-stats-json-size.md](N-homoplasy-stats-json-size.md)). With those, the current reader takes 0.4 s instead of 9.8 s for the `sc2/4500` statistics in the dev profile.

> [!IMPORTANT]
> **Decision required.** Options for the reader, alone or combined:
>
> - **Concrete reader**: return a concrete `BufReader<Decompressor>` from `open_file_or_stdin()` for both files and standard input, with standard input read through a `Decompressor` without compression, as `create_file_or_stdout()` already keeps `BufWriter<Compressor>` outermost. The boxed reader inside `Decompressor` is then called once per 256 KiB refill. 2 to 3 times faster, still streaming, no memory cost
> - **Parse from memory**: read the decompressed input into a buffer and parse it with `serde_json::from_slice()` behind `serde_stacker`. 6 to 16 times faster; peak memory grows by the size of the decompressed file
> - **Memory map**: map an uncompressed file and parse it with `from_slice()`. As fast as parsing from memory without a heap copy, but it needs `unsafe` and a new dependency (`memmap2`), the behavior is undefined if another process changes the file while it is mapped, and compressed files and standard input need another path
> - **Buffered stream reader**: use a serde_json build with a reader that scans a fixed buffer in bulk, such as serde-rs/json#1294. As fast as parsing from memory, with memory bounded by the buffer, but it needs a patched serde_json because the `Read` trait is sealed. The pull request is open without a maintainer review since 2025-10-29; an earlier change that buffered inside `from_reader()` (serde-rs/json#1007) was closed because a caller may read the stream after the JSON value
