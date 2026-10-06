# Tree JSON outputs are written without indentation

The JSON tree outputs (Auspice JSON, graph JSON, and UShER MAT JSON) are written without indentation or line breaks. The other JSON outputs, such as the GTR model, the clock model, the coalescent, the augur node data, and the `homoplasy` statistics, keep 2-space indentation.

**Type**: Output-format choice for outputs that have no v0 counterpart.

**Status**: Approved by the maintainer team.

**v1**: `fn write_graph_outputs()`, `fn write_tree_formats()`, and `fn write_mat_outputs()` in [`packages/app-output/src/tree_output.rs`](../../packages/app-output/src/tree_output.rs) write these files with `JsonPretty(false)`. Every other JSON output is written with `JsonPretty(true)`.

## Decision

Write the JSON files that describe a tree without whitespace, and keep the indentation of the other JSON outputs.

## Rationale

- In an Auspice tree, the node object and its `children` array each add a nesting level per tree level, so with indentation the size of a file grows with the number of nodes times the tree depth. Indented Auspice files are 91 to 97% whitespace: the Auspice JSON of `sc2/4500` is 47.8 MB with indentation and 1.3 MB without, and a `timetree.auspice.json` of `mpox/clade-ii/500` is 199 MB and 17.8 MB
- Every reader parses the whitespace byte by byte. With the stream reader of `json_read_file()` in the dev profile, the Auspice JSON of `sc2/4500` takes 3.7 s to read with indentation and 0.26 s without ([kb/issues/M-io-json-read-from-reader-slow.md](../issues/M-io-json-read-from-reader-slow.md)); writing the `mpox/clade-ii/500` tree takes 5.9 billion CPU instructions with indentation and 0.66 billion without
- People rarely read tree files directly; they open them in Auspice or another viewer. The smaller JSON outputs, such as model parameters, are read directly, so they stay indented
- A fixed rule per output kind keeps the format of each output independent of the dataset size and needs no extra pass to measure the size

## Limits

- An indented copy of a tree file needs a formatter, for example `jq . file.json`
- augur `export v2` uses a size rule instead: since augur 24.0.0 it writes without indentation when the indented JSON would exceed 5 MB, measured by serializing the data once into a byte counter (`json_size()` in `augur/io/json.py`). Small Auspice files of TreeTime are therefore compact where augur indents them; Auspice reads both forms
- JSON outputs that are not trees but grow with the dataset keep their indentation; the `homoplasy` statistics file is 38 to 52% whitespace on large datasets ([kb/issues/N-homoplasy-stats-json-size.md](../issues/N-homoplasy-stats-json-size.md))
