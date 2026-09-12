# treetime-io

I/O for phylogenetic data formats: trees, sequences, and metadata.

Generic serialization helpers (JSON, YAML) live in `treetime-utils::io`.

## Formats

| Module                | Format                                                       | Read | Write |
| --------------------- | ------------------------------------------------------------ | ---- | ----- |
| `nwk`                 | Newick tree format                                           | yes  | yes   |
| `nex`                 | Nexus tree format                                            | -    | yes   |
| `fasta`               | FASTA sequences                                              | yes  | yes   |
| `auspice`             | Auspice v2 JSON (Nextstrain visualization)                   | yes  | yes   |
| `phyloxml`            | PhyloXML tree format                                         | yes  | yes   |
| `usher_mat`           | UShER MAT protobuf                                           | yes  | yes   |
| `csv`                 | CSV/TSV/SSV with bounded content-based delimiter detection   | yes  | yes   |
| `dates_csv`           | Date metadata from CSV/TSV/SSV (year fractions, date ranges) | yes  | -     |
| `discrete_states_csv` | Discrete trait metadata from CSV/TSV/SSV                     | yes  | -     |
| `graphviz`            | Graphviz DOT format                                          | -    | yes   |

## Utilities

| Module            | Purpose                                                        |
| ----------------- | -------------------------------------------------------------- |
| `auspice_types`   | Auspice v2 JSON schema types                                   |
| `concat`          | Reader adaptor that concatenates multiple readers sequentially |
| `parse_delimited` | Split input by byte delimiter into string chunks               |

## Key types

- `FastaRecord` - sequence name, description, sequence data, index
- `DateConstraint` - parsed date with original input string, distinguishing exact, uncertain, and range inputs
- `CsvStructWriter` / `CsvStructFileWriter` - typed CSV writing via serde

## Tree I/O values

Readers return the parsed graph alongside per-node and per-edge value maps (`NwkParse`: names, confidence, branch lengths). Writers take those value maps plus optional `NodeCommentProvider`s for node annotations, so serialization reads each node's and edge's data from the threaded value maps.

## API patterns

File functions accept `impl AsRef<Path>` and support stdin/stdout via `treetime-utils` file utilities. Each format provides `_file`, `_str`, and raw reader/writer variants where applicable.

```rust
use treetime_io::fasta::{FastaRecord, read_fasta_file};
use treetime_utils::io::json::{json_read_file, json_write_file, JsonPretty};

// FASTA
let records: Vec<FastaRecord> = read_fasta_file("sequences.fasta", &alphabet)?;

// JSON (from treetime-utils)
let config: Config = json_read_file("config.json")?;
json_write_file("output.json", &config, JsonPretty(true))?;
```

## Notes

- File readers transparently decompress `.xz`, `.gz`, `.bz2`, and `.zst` files.
- Metadata readers select a delimiter from `--metadata-delimiters` by parsing a bounded sample of the decompressed header. The logical `.csv`, `.tsv`, or `.ssv` extension resolves ambiguous samples.
