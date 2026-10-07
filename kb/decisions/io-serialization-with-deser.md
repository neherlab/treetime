# Serialization uses deser

Serialized types of the workspace derive `deser::Serialize`, `deser::Deserialize`, or both, only for the directions that code writes or reads, and JSON, CSV, YAML, XML, and query strings are read and written with the deser format crates. serde stays only where a library requires `serde_json::Value`.

**Type**: Library choice for an implementation concern without a v0 counterpart.

**Status**: Proposed, awaiting approval by the maintainer team.

**v1**:

- JSON: `json_read()`, `json_write()`, and the other helpers of [`packages/treetime-utils/src/io/json.rs`](../../packages/treetime-utils/src/io/json.rs) use deser-json, with `deser_path::PathLayer` so that errors name the path of the offending value
- Tables: [`packages/treetime-io/src/csv.rs`](../../packages/treetime-io/src/csv.rs) uses deser-csv for typed rows, metadata records, and the delimiter probe
- YAML: `yaml_read_str()` and `yaml_value_read_str()` in [`packages/app-commands/src/yaml.rs`](../../packages/app-commands/src/yaml.rs) use deser-yaml for configuration and settings files; writing YAML keeps the `saphyr` emitter
- Server: request and response bodies, query strings (deser-urlencoded), and path parameters of `packages/app-server` use deser
- Dynamic JSON: schemars, jsonschema, aide, and minijinja work on `serde_json::Value`. `JsonValue` and `SparseConfig` in [`packages/app-commands/src/json_value.rs`](../../packages/app-commands/src/json_value.rs) and the `As<Value, Serde>` adapter of deser-serde carry these values across

## Decision

Replace serde derives and the serde format crates (serde_json for typed values, serde-saphyr, csv, serde_with, serde_stacker) with deser 0.9.0 and its format crates. Keep `serde_json::Value` for the documents that schema and template libraries take, bridged with deser-serde.

## Rationale

Release profile, medians of repeated runs, files in the page cache, on a thread with a 2 MiB stack; "before" is the same workspace with serde:

| Operation                                                 | Before  | deser   |
| --------------------------------------------------------- | ------- | ------- |
| Read `sc2/4500` Auspice JSON, indented (47.8 MB)          | 0.251 s | 0.051 s |
| Read `mpox/clade-ii/500` Auspice JSON, indented (199 MB)  | 1.08 s  | 0.287 s |
| Read `mpox/clade-ii/500` Auspice JSON, compact (17.8 MB)  | 0.210 s | 0.175 s |
| Read `sc2/4500` `homoplasy` statistics, indented (125 MB) | 0.703 s | 0.320 s |
| Write `mpox/clade-ii/500` Auspice JSON, indented          | 0.720 s | 0.246 s |
| Write `mpox/clade-ii/500` Auspice JSON, compact           | 0.068 s | 0.053 s |

- **Reading streams at the speed of parsing from memory.** The stream reader of deser-json scans its input buffer in bulk, so `json_read_file()` takes as long as parsing the whole file from a string, with memory bounded by the parsed value: the `sc2/4500` statistics file reads in 0.317 s from the stream (14 MB peak) and 0.301 s from memory (133 MB peak). The serde_json stream reader handled one byte per call and was 3 to 7 times slower than its slice parser
- **Fast in the dev profile.** The dev server and the tests use the dev profile, where the serde_json stream reader was 4 to 6 times slower than in the release profile. With deser, the `mpox/clade-ii/500` Auspice JSON reads in 0.25 s instead of 5.9 s, and the `sc2/4500` statistics file in 0.28 s instead of 3.7 s
- **No stack limit on nesting.** The deser drivers keep the nesting on the heap. A 50,000-level Auspice tree reads and writes on a 2 MiB stack; with serde, writing it overflowed the stack, and reading it needed `serde_stacker` and 670 MB of peak memory instead of 199 MB
- **Errors name the setting.** Errors of JSON, YAML, and query-string input name the path of the offending value, for example `(path: paths)`
- **Same outputs.** The generated JSON schemas, the OpenAPI document, the TypeScript client, and the CLI reference are unchanged. The full smoke matrix, compared with the same workspace before the change, gives identical outputs for every case that both versions complete; the augur node data of config-file runs differs only in the absolute input paths, which contain the folder of each snapshot

## Limits

- **Rebuilds after an edit are slower.** The deser derive generates more code per type, so a rebuild that reaches `app-commands` takes 4 to 5 s longer ([kb/issues/M-build-deser-derives-slow-down-rebuilds.md](../issues/M-build-deser-derives-slow-down-rebuilds.md)). A clean release build takes about as long as before
- **No flattened untagged or adjacently tagged enums.** deser cannot flatten them, so the run record, the run event, and the Auspice CDS types use internally tagged enums and optional fields that keep the same JSON
- **serde stays in the dependency tree.** schemars, jsonschema, aide, axum, napi, and minijinja depend on serde, and the schema transforms take defaults and skipped optional fields from deser because schemars reads defaults only through serde `Serialize`
- **YAML errors have no code excerpt.** Settings and configuration errors give the line, the column, and the setting path, without the source excerpt that serde-saphyr printed
- **YAML 1.2 with YAML 1.1 boolean words.** deser-yaml reads YAML 1.2; a layer keeps plain `yes`, `no`, `on`, `off`, `y`, `n` as booleans. The YAML 1.1 mode of deser-yaml is not used, because its float syntax needs a decimal point and reads `1e-12` as text
- **Young dependency.** deser 0.9.0 has one maintainer and was released in 2026-09; the format crates follow its version
