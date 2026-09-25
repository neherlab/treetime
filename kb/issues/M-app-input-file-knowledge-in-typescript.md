# The app decides input-file names, formats and check requests in TypeScript

The TypeScript apps are a user interface around the Rust core (see "Apps" in `docs/dev/developer_guide.md`). Settings, command lines, YAML and pre-run checks come from Rust, but three pieces of input-file knowledge remain in `packages/app-ui`:

- **Check-inputs request**: `inputFactsRequest` in `packages/app-ui/src/settings/inputs.ts` maps a draft configuration to a `CheckInputsRequest` by setting name (`tree`, `metadata`, `alignment`, `metadata_id_columns`, `metadata_delimiters`, `date_column`). The config types know these settings; a Rust operation that takes the command and its configuration would read them the way the command does
- **Dataset file names**: `DATASET_FILES` in the same file names the files of an example dataset (`tree.nwk`, `aln.fasta.xz`, `aln.fasta`, `metadata.tsv`, `metadata.csv`). `packages/app-datasets/src/lib.rs` discovers datasets by `tree.nwk` on its own, so the two lists can diverge
- **Accepted extensions**: `INPUT_SLOT_INFO` in `packages/app-ui/src/settings/commands.ts` lists the file extensions the file picker offers per input. The readers in `treetime-io` decide which formats and compressions are read (`guess_compression_from_filepath`), so a new format or compression needs a second edit in TypeScript

## Impact

- A new metadata setting, dataset file name, or compression format works in the CLI and is missed by the app until the TypeScript copy is updated

## Potential solutions

- O1. Let `check-inputs` take the command and its configuration, and read the input settings with the config types
- O2. Return the input files of each dataset from `app-datasets` with the input kind they fill
- O3. Add the accepted extensions of each input kind to the setting catalog, from the readers of `treetime-io`
