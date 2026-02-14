# treetime-io

I/O for phylogenetic data formats: trees, sequences, metadata, and serialization.

## Formats

| Module                | Format                                                          | Read | Write |
| --------------------- | --------------------------------------------------------------- | ---- | ----- |
| `nwk`                 | Newick tree format                                              | yes  | yes   |
| `nex`                 | Nexus tree format                                               | -    | yes   |
| `fasta`               | FASTA sequences (with transparent xz/gz/bz2/zstd decompression) | yes  | yes   |
| `auspice`             | Auspice v2 JSON (Nextstrain visualization)                      | yes  | yes   |
| `phyloxml`            | PhyloXML tree format                                            | yes  | yes   |
| `usher_mat`           | UShER MAT protobuf                                              | yes  | yes   |
| `json`                | JSON via `serde_json` (with `serde_stacker` for deep nesting)   | yes  | yes   |
| `yaml`                | YAML via `serde_yaml`                                           | yes  | yes   |
| `csv`                 | CSV/TSV with auto-delimiter detection                           | yes  | yes   |
| `dates_csv`           | Date metadata from CSV/TSV (year fractions, date ranges)        | yes  | -     |
| `discrete_states_csv` | Discrete trait metadata from CSV/TSV                            | yes  | -     |
| `graphviz`            | Graphviz DOT format                                             | -    | yes   |

## Utilities

| Module            | Purpose                                                        |
| ----------------- | -------------------------------------------------------------- |
| `auspice_types`   | Auspice v2 JSON schema types                                   |
| `concat`          | Reader adaptor that concatenates multiple readers sequentially |
| `parse_delimited` | Split input by byte delimiter into string chunks               |

## Key types

- `FastaRecord` - sequence name, description, sequence data, index
- `DateOrRange` - year fraction or year fraction range parsed from date metadata
- `CsvStructWriter` / `CsvStructFileWriter` - typed CSV writing via serde
- `JsonPretty(bool)` - controls JSON pretty-printing

## Tree I/O traits

Tree readers/writers use trait bounds on graph node and edge types:

- `NodeFromNwk` / `EdgeFromNwk` - construct nodes and edges from Newick data
- `NodeToNwk` / `EdgeToNwk` - serialize nodes and edges to Newick
- `NodeToGraphviz` / `EdgeToGraphviz` - serialize to Graphviz DOT
- `AuspiceRead` / `AuspiceWrite` - convert to/from Auspice v2 JSON
- `UsherRead` - convert from UShER MAT protobuf
- `PhyloxmlToGraph` / `PhyloxmlDataToGraphData` - convert to/from PhyloXML

## API patterns

File functions accept `impl AsRef<Path>` and support stdin/stdout via `treetime-utils` file utilities. Each format provides `_file`, `_str`, and raw reader/writer variants where applicable.

```rust
use treetime_io::json::{json_read_file, json_write_file, JsonPretty};
use treetime_io::yaml::{yaml_read_file, yaml_write_file};
use treetime_io::fasta::{FastaRecord, read_fasta_file};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
struct Config {
    name: String,
    value: i32,
}

// JSON
let config: Config = json_read_file("config.json")?;
json_write_file("output.json", &config, JsonPretty(true))?;

// YAML
let config: Config = yaml_read_file(&Some("config.yaml"))?;
yaml_write_file("output.yaml", &config)?;
```

## Notes

- JSON reading disables recursion limits to handle deeply nested structures
- `json_or_yaml_write_file` determines format by file extension (`.yaml`/`.yml` for YAML, otherwise JSON)
- FASTA reading transparently decompresses `.xz`, `.gz`, `.bz2`, `.zst` files
- CSV delimiter is auto-detected from file extension (`.tsv` uses tab, `.csv` uses comma)
