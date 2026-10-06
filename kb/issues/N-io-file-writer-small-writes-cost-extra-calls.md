# Small writes through FileWriter cost extra calls

`FileWriter` in `packages/treetime-utils/src/io/file.rs` implements only `write()` and `flush()` of `std::io::Write`. A caller that writes many small pieces with `write_all()` goes through the default `write_all()` loop and a non-inlined call to `FileWriter::write()` and then `BufWriter::write()` for each piece. The pretty printer of serde_json writes the indentation of every line as one 2-byte `write_all()` per nesting level, so indented JSON outputs pay this cost for every level of every line.

## Evidence

Writing an indented `timetree.auspice.json` of `mpox/clade-ii/500` (199 MB) with `json_write_file()` takes 5.9 billion CPU instructions in the `profiling` profile. The same `serde_json::to_writer_pretty()` into a plain `BufWriter<File>` takes 3.4 billion instructions.

Tree JSON outputs are written without indentation ([kb/decisions/io-tree-json-outputs-without-indentation.md](../decisions/io-tree-json-outputs-without-indentation.md)), so the cost applies to the indented JSON outputs, such as the `homoplasy` statistics file, and to any other writer that makes many small writes.

## Fix direction

Forward `write_all()` of `FileWriter` to `BufWriter::write_all()`, which copies small pieces into its buffer without a loop.

> [!IMPORTANT]
> **Investigation required.** Measure the change on the indented statistics file of `sc2/4500` and on a compact Auspice file before and after.
