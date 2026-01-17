# treetime-io

I/O utilities for reading and writing JSON and YAML files with serde serialization.

## Supported Formats

- **JSON** - via `serde_json` with support for deeply nested structures (uses `serde_stacker` to avoid recursion limits)
- **YAML** - via `serde_yaml`

## API

### JSON

| Function                             | Description                                                             |
| ------------------------------------ | ----------------------------------------------------------------------- |
| `json_read_file(path)`               | Read and deserialize JSON from file                                     |
| `json_read_str(s)`                   | Deserialize JSON from string                                            |
| `json_read(reader)`                  | Deserialize JSON from any `std::io::Read`                               |
| `json_write_file(path, obj, pretty)` | Serialize and write JSON to file                                        |
| `json_write_str(obj, pretty)`        | Serialize to JSON string                                                |
| `json_write(writer, obj, pretty)`    | Serialize JSON to any `std::io::Write`                                  |
| `json_or_yaml_write_file(path, obj)` | Write JSON or YAML based on file extension                              |
| `is_json_value_null(t)`              | Check if value serializes to null (for `#[serde(skip_serializing_if)]`) |

### YAML

| Function                     | Description                               |
| ---------------------------- | ----------------------------------------- |
| `yaml_read_file(path)`       | Read and deserialize YAML from file       |
| `yaml_read_str(s)`           | Deserialize YAML from string              |
| `yaml_read(reader)`          | Deserialize YAML from any `std::io::Read` |
| `yaml_write_file(path, obj)` | Serialize and write YAML to file          |
| `yaml_write_str(obj)`        | Serialize to YAML string                  |
| `yaml_write(writer, obj)`    | Serialize YAML to any `std::io::Write`    |

## Usage

```rust
use treetime_io::json::{json_read_file, json_write_file, JsonPretty};
use treetime_io::yaml::{yaml_read_file, yaml_write_file};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
struct Config {
    name: String,
    value: i32,
}

// Read JSON
let config: Config = json_read_file("config.json")?;

// Write JSON (pretty-printed)
json_write_file("output.json", &config, JsonPretty(true))?;

// Read YAML
let config: Config = yaml_read_file(&Some("config.yaml"))?;

// Write YAML
yaml_write_file("output.yaml", &config)?;
```

## Notes

- File functions support stdin/stdout via `treetime-utils` file utilities
- JSON reading disables recursion limits to handle deeply nested structures
- `json_or_yaml_write_file` determines format by extension (`.yaml`/`.yml` for YAML, otherwise JSON)
