# Implement three-tier output selection

Replace the hardcoded four-file output with a three-tier selection system: bulk directory, output selection filter, per-file path overrides.

## Scope

- Define `OutputSelection` enum covering tree formats (nwk, nwk-annotated, nwk-nhx, nexus, nexus-annotated, nexus-nhx, auspice, phyloxml, mat-pb, mat-json, graph-json, dot) and non-tree outputs (augur-node-data, fasta, clock-model, traits-csv, gtr, confidence)
- Replace `OutputArgs.outdir` with three tiers:
  - `--output-all` (`-O`): bulk directory with command-specific default basename
  - `--output-selection` (`-s`): comma-separated filter for which outputs `--output-all` produces
  - Per-file `--output-tree-<format>=<path>` and `--output-<non-tree>=<path>` overrides
- Resolution logic: `--output-all` fills defaults via `get_or_insert`, `--output-selection` restricts, per-file flags always win, explicit per-file outside selection still honored
- Replace `write_graph_files_with_options` at `graph.rs:18-53` with format-selection-aware writer
- Wire into all 6 commands (ancestral, timetree, optimize, mugration, clock, prune)
- Validate requested formats against per-command adapter trait availability (error on unavailable format)
- Default tree formats: `nwk,nexus` (intentionally drops graph-json and dot from current defaults)

## Locations

- OutputArgs: `packages/treetime/src/commands/shared/output.rs`
- Shared writer: `packages/treetime-io/src/graph.rs:18-53`
- NwkStyle enum from ticket #1: used to distinguish nwk/nwk-annotated/nwk-nhx variants
- Command call sites: 6 `run.rs` files

## Tests

### Happy paths -- tier 1 (bulk directory)

- `--output-all dir/` creates directory, produces `dir/stem.nwk` and `dir/stem.nexus` (defaults)
- Default basename per command: `timetree` for timetree, `annotated_tree` for ancestral/optimize/mugration, `rerooted` for clock, `pruned_tree` for prune
- Creates parent directories if missing

### Happy paths -- tier 2 (selection filter)

- `--output-all dir/ --output-selection=nwk` produces only `.nwk`
- `--output-all dir/ --output-selection=nwk,nwk-annotated,nexus` produces three files
- `--output-all dir/ --output-selection=all` produces every format available for that command
- `--output-selection` without `--output-all` is ignored (Nextclade behavior: `requires = "output_all"`)

### Happy paths -- tier 3 (per-file overrides)

- `--output-tree-nwk=out.nwk` produces only that file at that path, no `--output-all` needed
- `--output-tree-nexus=out.nexus --output-tree-nwk=out.nwk` produces both
- `--output-all dir/ --output-tree-nwk=custom.nwk`: NWK at custom path, nexus at default path (per-file wins)
- `--output-all dir/ --output-selection=nwk --output-tree-nexus=extra.nexus`: selection says nwk only, but explicit per-file nexus also produced
- Non-tree: `--output-augur-node-data=node.json` produces augur JSON (ancestral, timetree)
- `-` as path: writes to stdout (uncompressed)

### Happy paths -- per-command

- Each of 6 commands produces correct output with `--output-all dir/`
- Non-tree outputs included: augur node data (ancestral, timetree), FASTA (ancestral), clock model (timetree, clock), traits CSV (mugration), GTR JSON

### Happy paths -- compression

- `.gz` extension: gzip-compressed output
- `.xz` extension: xz-compressed output
- `.zst` extension: zstd-compressed output
- `.bz2` extension: bzip2-compressed output

### Edge cases

- No output flags at all -> error (at least one required)
- `--output-all dir/ --output-selection=` (empty selection) -> error
- `--output-selection=auspice` on clock command -> error: format not available (no AuspiceWrite impl for NodeClock)
- `--output-selection=mat-pb` on mugration -> error: no UsherWrite impl, no mutation data
- Path with spaces: `--output-tree-nwk="my dir/tree.nwk"` -> works
- Relative path: `--output-tree-nwk=../out.nwk` -> resolves correctly
- Same path for two different flags -> last write wins or error? Document behavior
- `--output-all` to existing directory with existing files -> overwrites without warning (standard behavior)

### Pathological

- Output to read-only directory -> IO error propagated with context
- Output path with 260+ chars (Windows limit) -> OS error propagated
- `--output-selection` with unknown format name -> error listing valid options
- `--output-selection` with duplicates: `nwk,nwk` -> deduplicated or error? Document

### Regression

- Default `--output-all dir/` produces same tree content as current `write_graph_files_with_options` (topology, branch lengths, names identical in nwk and nexus)
- graph-json and dot no longer produced by default (intentional change, documented)

## Related issues

- Source: [kb/issues/N-io-write-graph-files-missing-formats.md](../issues/N-io-write-graph-files-missing-formats.md)
- Source: [kb/issues/N-io-nexus-writer-no-annotation-advantage.md](../issues/N-io-nexus-writer-no-annotation-advantage.md)
- Proposal: [kb/proposals/output-format-selection.md](../proposals/output-format-selection.md)
