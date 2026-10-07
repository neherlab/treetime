# Small writes through FileWriter cost extra calls

`FileWriter` in `packages/treetime-utils/src/io/file.rs` implements only `write()` and `flush()` of `std::io::Write`. A caller that writes many small pieces with `write_all()` goes through the default `write_all()` loop and a non-inlined call to `FileWriter::write()` and then `BufWriter::write()` for each piece.

The JSON and CSV outputs do not pay this cost: deser-json and deser-csv serialize into their own buffer and pass it to the writer in pieces of at least 8 KiB (`DEFAULT_BUFFER_LIMIT` of `deser_core::io::Writer`). The Newick and Nexus writers in `packages/util-newick/src/write/newick.rs` and `packages/util-newick/src/nexus/write.rs` do not pay it either: they format each tree into a `String` and pass it to the writer with one `write_all()`. The writers that format their output piece by piece with `write!()` and `write_all()` still pay it, such as the FASTA writer in `packages/treetime-io/src/fasta.rs` and the Graphviz writer in `packages/treetime-io/src/graphviz.rs`.

## Fix direction

Forward `write_all()` of `FileWriter` to `BufWriter::write_all()`, which copies small pieces into its buffer without a loop.

> [!IMPORTANT]
> **Investigation required.** Measure the CPU instructions of writing a large FASTA file and a large Graphviz file through `FileWriter` before and after the change.
