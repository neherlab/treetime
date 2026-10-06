# Indented JSON output makes tree files 10 to 37 times larger

Every JSON output file is written with `JsonPretty(true)`, which uses the pretty printer of serde_json with 2 spaces per nesting level (for example `TreeWriteKind::Auspice` in `packages/app-output/src/tree_output.rs` and the statistics file in `packages/app-commands/src/commands/homoplasy/run.rs`). In an Auspice tree, the node object and its `children` array each add a nesting level per tree level, so the indentation of a line grows with the depth of its node, and the file size grows with the number of nodes times the tree depth. Large Auspice files are mostly spaces, and every program that reads them, including the results pages of the apps, parses those spaces byte by byte.

## Evidence

| File                                           | Indented | Whitespace | Without whitespace |
| ---------------------------------------------- | -------- | ---------- | ------------------ |
| Auspice JSON of `sc2/4500` (tree depth 157)    | 47.8 MB  | 97.4%      | 1.3 MB             |
| `timetree.auspice.json` of `mpox/clade-ii/500` | 198.9 MB | 91.0%      | 17.8 MB            |
| `timetree.auspice.json` of `rsv/a/2000`        | 19.4 MB  | 94.8%      | 1.0 MB             |
| `homoplasy.stats.json` of `sc2/4500`           | 125.0 MB | 38.2%      | 77.3 MB            |
| `homoplasy.stats.json` of `mpox/clade-ii/2000` | 226.7 MB | 52.4%      | 107.9 MB           |

The mean indentation in the Auspice JSON of `sc2/4500` is 332 spaces per line.

- **Reading**: in the dev profile, `json_read_file()` takes 3.7 s for the indented Auspice JSON of `sc2/4500` and 0.26 s for the same tree without whitespace, and 15.7 s and 1.83 s for `mpox/clade-ii/500` ([M-io-json-read-from-reader-slow.md](M-io-json-read-from-reader-slow.md))
- **Writing**: `json_write_file()` writes the `mpox/clade-ii/500` Auspice tree with 5.9 billion CPU instructions indented and 0.66 billion without indentation (about 0.6 s and 0.06 s in the `profiling` profile). The pretty printer writes the indentation as one 2-byte `write_all()` per nesting level. `FileWriter` in `packages/treetime-utils/src/io/file.rs` implements only `write()`, so each piece goes through the default `write_all()` loop and a call to `BufWriter::write()`; a plain `BufWriter<File>` writes the same indented file with 3.4 billion instructions
- **Disk and transfer**: the run folders of the apps and the downloads of their outputs carry the whitespace

augur `export v2` writes without indentation when the indented JSON would be larger than 5 MB, since augur 24.0.0 (nextstrain/augur#1352); `--minify-json` and `--no-minify-json` override the threshold. augur does not know the size in advance: `json_size()` in `augur/io/json.py` serializes the data once with indentation into a stream that only counts bytes, and `write_json()` then serializes it again into the file. The threshold applies to the indented size, so a deep tree that is small without indentation is still written without indentation. Auspice reads both forms.

The same counting pass in Rust, `serde_json::to_writer()` into a `Write` that only adds up the buffer lengths, costs 0.14 billion instructions for the `mpox/clade-ii/500` Auspice tree in the `profiling` profile, the same for the indented and the compact form, against 0.66 billion for writing the compact file and 5.9 billion for writing the indented file.

A change of the format changes every affected JSON output byte for byte, so the smoke baseline and the JSON reference fixtures change with it.

> [!IMPORTANT]
> **Decision required.** Options for the JSON output files:
>
> - **No indentation for all JSON outputs**: smallest files and fastest reads and writes; the files are harder to read without a formatter such as `jq`
> - **No indentation above a size threshold**: the augur rule (indented size above 5 MB), so small outputs stay readable; it needs a counting pass before each write, and the format of an output then depends on the size of the dataset
> - **No indentation for tree files only**: Auspice JSON and the other files whose size grows with the tree depth; the statistics file keeps its indentation and stays about 2 times larger than needed
> - **Keep the indentation**: the costs above remain; a narrower indentation unit only shrinks the files by a constant factor
>
> Independent of the format: `FileWriter` can forward `write_all()` to its `BufWriter`, which removes the extra calls for every small write.
