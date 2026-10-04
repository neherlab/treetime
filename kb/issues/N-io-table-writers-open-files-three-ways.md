# Table writers open their output files three ways

Every CSV and TSV output goes through `CsvStructWriter` in `packages/treetime-io/src/csv.rs`, which quotes fields that contain the delimiter, a quote, or a line break. The writers differ in how they open the destination and attach error context:

- `CsvStructFileWriter` (`.clock.csv` in `packages/app-output/src/rtt.rs`, `.tracelog.csv` in `packages/app-output/src/timetree_trace.rs`): opens through `create_file_or_stdout()`; write and flush errors carry no file path
- `write_file_with()` + `CsvStructWriter` (`.confidence.tsv` in `packages/app-output/src/confidence.rs`, `.traits.csv` and `.confidence.csv` in `packages/app-output/src/mugration_result.rs`): write and flush errors carry the context `When writing file '<path>'`
- `create_file_or_stdout()` + `CsvStructWriter` + `FileWriter::finish()` (`.coalescent.tsv` and `.coalescent.csv` in `packages/app-output/src/coalescent.rs`): write and flush errors carry no file path

The escaping is the same in all three, but a write failure reports the failed path for some tables and not for others, and each new table writer copies one of the three patterns.

## Constraints

- The timetree trace writes one record per iteration while inference runs, so it needs a writer that stays open across calls. The other tables are written once from complete data.

## Decision axes

- **One entry point**: make `write_file_with()` + `CsvStructWriter` the only pattern for complete tables and keep `CsvStructFileWriter` for streaming, or give `CsvStructFileWriter` the same path context and use it everywhere
- **Ownership**: keep table writing in `app-output` per format, or move it to the output dispatch when the tree-output data model is redesigned
