# Write augur node data JSON from mugration command

Produce augur-compatible node data JSON output from `run_mugration()`. This is the treetime equivalent of `augur traits` output. Written to `<outdir>/mugration.augur-node-data.json` by default.

## CLI argument

```rust
/// Write augur-compatible node data JSON to this path
///
/// Contains per-node discrete trait assignments, confidence profiles, entropy,
/// the inferred substitution model, and branch state-change labels. The output
/// is compatible with augur export v2 --node-data for Nextstrain pipeline
/// integration.
#[cfg_attr(feature = "clap", clap(long))]
#[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
#[cfg_attr(feature = "clap", clap(default_value = "<outdir>/mugration.augur-node-data.json"))]
pub output_augur_node_data: PathBuf,
```

## Output fields

Top-level: `generated_by`, `models.<column>` (rate, alphabet including "?" missing marker, equilibrium_probabilities excluding missing, transition_matrix excluding missing), optional `branches` (branch labels at state-change nodes).

Per-node: `<column>` (inferred state string), `<column>_confidence` (state-to-probability map, sorted descending, filtered > 0.001), `<column>_entropy` (Shannon entropy: `-sum(p * ln(p + 1e-12))`).

## Data availability

`MugrationResult` carries the partition. All needed data is accessible: `traits.assignments` (trait values), `confidence.rows[].profile` (posteriors), `partition.states` (state names), `partition.gtr()` (model parameters). No partition scope problem.

Branch labels require comparing parent/child trait assignments. Arrow character is Unicode U+2192.

## Alphabet sizing asymmetry

`models.<column>.alphabet` includes the missing data marker at the end (n_states+1 elements). `equilibrium_probabilities` and `transition_matrix` exclude it (n_states elements). This is intentional.

## Dependencies

- `io-implement-augur-node-data-json-types.md` (serde types)

## Reference

- Format specification: [../reports/augur-node-data-json.md](../reports/augur-node-data-json.md) section "`augur traits` node data JSON"
- Augur transformation code: [`traits.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/traits.py)
