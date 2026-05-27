# Write augur node data JSON from ancestral command

Produce augur-compatible node data JSON output from `run_ancestral_reconstruction()`. This is the treetime equivalent of `augur ancestral` output. Written to `<outdir>/ancestral.augur-node-data.json` by default.

## CLI argument

```rust
/// Write augur-compatible node data JSON to this path
///
/// Contains per-node nucleotide mutations, reconstructed sequences, alignment
/// mask, genome annotations, and the reference sequence. The output is
/// compatible with augur export v2 --node-data for Nextstrain pipeline
/// integration.
#[cfg_attr(feature = "clap", clap(long))]
#[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
#[cfg_attr(feature = "clap", clap(default_value = "<outdir>/ancestral.augur-node-data.json"))]
pub output_augur_node_data: PathBuf,
```

## Output fields

Top-level: `generated_by`, `annotations.nuc` (start=1, end=alignment_length, strand="+", type="source"), `reference.nuc` (root sequence), `mask` (per-position ambiguity string).

Per-node: `muts` (mutation strings from `edge_subs()`, 1-based via `Sub::Display`), `sequence` (full reconstructed sequence).

AA fields (`aa_muts`, `aa_sequences`, gene annotations) are out of scope. See [../proposals/node-data-json-aa-reconstruction.md](../proposals/node-data-json-aa-reconstruction.md).

## Partition scope

The partition holds the needed data (mutations, sequences, root sequence, alignment length) but is local to `run_ancestral_reconstruction()`. Two options:

- Write inside the run function scope before partition drops (alongside existing FASTA/Newick output)
- Extend `AncestralResult` to carry the needed data

The choice affects separation of concerns. The report discusses tradeoffs in the "Partition access constraint" subsection.

## Dependencies

- `io-implement-augur-node-data-json-types.md` (serde types)
- `ancestral-implement-mask-computation.md` (mask field + mutation filtering)

## Reference

- Format specification: [../reports/augur-node-data-json.md](../reports/augur-node-data-json.md) section "`augur ancestral` node data JSON"
- Augur transformation code: [`ancestral.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/ancestral.py)
