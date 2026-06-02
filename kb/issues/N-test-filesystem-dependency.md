# Tests with unnecessary filesystem dependency

Several unit tests use `tempfile`/`tempdir` or read `data/` files when the underlying logic could be tested in-memory. This adds filesystem coupling, slows test execution, and makes tests fragile to working-directory or file-availability changes.

12 test files across 3 categories have filesystem dependency. After excluding golden master tests (exempt), smoke/integration tests (inherently FS-dependent), and tests of filename logic, 25 individual tests across 7 files are convertible to in-memory.

## Convertible: augur node data serialization (4 files, 17 tests)

All four `test_augur_node_data.rs` files follow the same pattern: `NamedTempFile::new()` then `write_augur_node_data_json(..., tmp.path())` then `std::fs::read_to_string(tmp.path())`. The writers internally call a `build_*` function to produce a struct, then serialize via `json_write_file`. Tests only care about the JSON string.

Conversion: call `build_*` directly, serialize with `json_write_str(&data, JsonPretty(true))`.

- `commands/ancestral/__tests__/test_augur_node_data.rs` - 4 tests via `helpers::write_json`. `pub fn build_augur_node_data_json` already exists at `commands/ancestral/augur_node_data.rs:35:`.
- `commands/mugration/__tests__/test_augur_node_data.rs` - 2 tests via `run_and_serialize`. Build logic is inlined in the writer body; extract a `build_augur_node_data_json` function from `commands/mugration/augur_node_data.rs`.
- `commands/optimize/__tests__/test_augur_node_data.rs` - 5 unit tests via `helpers::write_and_read` / `helpers::write_json`. Extract a `build_augur_node_data_json` from `commands/optimize/augur_node_data.rs`.
- `commands/timetree/output/__tests__/test_augur_node_data.rs` - 6 tests via `helpers::SampleCase::write_json`. The writer's core is `json_write(W: Write, ...)` so it already accepts a generic writer; tests could pass a `Vec<u8>` buffer. Or call the build function + `json_write_str`.

## Convertible: auspice JSON (1 file, 8 tests)

`commands/timetree/output/__tests__/test_auspice.rs` - 8 of 9 tests create `tempdir()`, write auspice JSON, read it back. `pub fn build_timetree_auspice` already exists at `timetree/output/auspice.rs:20:` and returns `Result<AuspiceTree>`. Tests can call it directly and assert on the struct. The two error-rejection tests (`test_auspice_rejects_nan_div`, `test_auspice_rejects_infinite_time`) already fail inside `build_timetree_auspice` before any file write.

Exception: `test_auspice_output_file_is_valid_json` genuinely tests file creation and should keep the tempdir.

## Convertible: data/ read-only (2 files, 3 tests, partial)

- `ancestral/__tests__/test_python_parity.rs` - 1 test reads `data/flu/h3n2/20/{tree.nwk,aln.fasta.xz}`. The tree is 1.5 KB (one line) and the alignment is 20 short sequences. Both can be embedded as `const` string literals, parsed with `nwk_read_str` / `read_many_fasta_str`. The other 6 tests in the file already use inline constants.
- `clock/__tests__/test_clock_dengue100.rs` - 2 tests read `data/dengue/100/{tree.nwk,metadata.tsv}`. The tree is 3.7 KB and the metadata is 101 lines (~2 KB). Embeddable as const, but the test's value is tied to this specific biological dataset (force_positive_rate behavior with all 198 root positions negative). Conversion changes character from "real dataset regression" to "inline fixture regression." Lower priority.

## Not convertible (must remain FS-dependent)

### Smoke/integration tests (3 files, 7 tests)

Full command e2e tests that exercise the entire pipeline including file I/O. Inherently FS-dependent.

- `commands/ancestral/__tests__/test_smoke_gtr_iterations.rs` - 2 tests: `run_ancestral_reconstruction` on `data/flu/h3n2/20/`, writes to `tmp/`, asserts output file existence and `result.gtr` properties.
- `commands/ancestral/__tests__/test_smoke_sample_from_profile.rs` - 2 tests: sampling reproducibility via file content comparison, output file existence. (1 additional rejection test is already essentially in-memory.)
- `commands/timetree/__tests__/test_pipeline.rs` - 1 test: full timetree pipeline, reads back tracelog CSV for convergence metrics, checks output file existence and size.

### E2e tests embedded in augur_node_data files (2 files, 4 tests)

- `commands/ancestral/__tests__/test_augur_node_data.rs` - 3 `reconstruct_json` tests: write synthetic tree+alignment to tempdir, run full `run_ancestral_reconstruction`, read back output JSON.
- `commands/optimize/__tests__/test_augur_node_data.rs` - 1 `end_to_end` test: reads `data/flu/h3n2/20/`, runs full optimize command, reads back output JSON.

### Filename logic tests (1 file, 2 tests)

- `gtr/__tests__/test_write_gtr_json.rs` - 2 tests: verify filename construction (qualified vs unqualified, non-overwriting). Inherently path-based.

### Large dataset convergence (1 file, 1 test)

- `optimize/__tests__/test_convergence_sc2.rs` - `test_convergence_sc2_sparse_converges_on_sc2_2844`: reads 149 KB tree + ~85 MB alignment. Not embeddable. The convergence bug is triggered by specific dataset properties. (The second test in this file uses `data/flu/h3n2/20/` which is small enough to embed but shares the fixture with `test_python_parity`.)

## Locations

All paths relative to `packages/treetime/src/`.

| File                                                             | Tests                             | Category        |
| ---------------------------------------------------------------- | --------------------------------- | --------------- |
| `commands/ancestral/__tests__/test_augur_node_data.rs`           | 4 convertible, 3 e2e              | augur node data |
| `commands/mugration/__tests__/test_augur_node_data.rs`           | 2 convertible                     | augur node data |
| `commands/optimize/__tests__/test_augur_node_data.rs`            | 5 convertible, 1 e2e              | augur node data |
| `commands/timetree/output/__tests__/test_augur_node_data.rs`     | 6 convertible                     | augur node data |
| `commands/timetree/output/__tests__/test_auspice.rs`             | 8 convertible, 1 keep             | auspice         |
| `ancestral/__tests__/test_python_parity.rs`                      | 1 convertible                     | data/ read      |
| `clock/__tests__/test_clock_dengue100.rs`                        | 2 convertible (lower priority)    | data/ read      |
| `commands/ancestral/__tests__/test_smoke_gtr_iterations.rs`      | 2 keep                            | smoke           |
| `commands/ancestral/__tests__/test_smoke_sample_from_profile.rs` | 2 keep                            | smoke           |
| `commands/timetree/__tests__/test_pipeline.rs`                   | 1 keep                            | smoke           |
| `gtr/__tests__/test_write_gtr_json.rs`                           | 2 keep                            | filename logic  |
| `optimize/__tests__/test_convergence_sc2.rs`                     | 1 keep (sc2), 1 convertible (flu) | data/ read      |
