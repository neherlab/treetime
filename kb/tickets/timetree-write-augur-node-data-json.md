# Write augur node data JSON from timetree command

Produce augur-compatible node data JSON output from `run_timetree_estimation()`. This is the treetime equivalent of `augur refine` output. Written to `<outdir>/timetree.augur-node-data.json` by default.

## CLI argument

```rust
/// Write augur-compatible node data JSON to this path
///
/// Contains per-node dates, branch lengths, clock model parameters, confidence
/// intervals, and divergence metrics. The output is compatible with augur
/// export v2 --node-data for Nextstrain pipeline integration.
#[cfg_attr(feature = "clap", clap(long))]
#[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
#[cfg_attr(feature = "clap", clap(default_value = "<outdir>/timetree.augur-node-data.json"))]
pub output_augur_node_data: PathBuf,
```

## Output fields

Top-level: `generated_by`, `clock` (rate, intercept, rtt_Tmrca, optional cov and rate_std), `alignment` (input path or null), `input_tree` (input path).

Per-node: `branch_length`, `confidence` (input tree bootstrap), `numdate`, `clock_length`, `mutation_length`, `raw_date` (original string from metadata), `date` (resolved YYYY-MM-DD), `date_inferred` (bool), optional `num_date_confidence` ([lower, upper]).

## Partition scope

Same constraint as ancestral: partition is local to `run_timetree_estimation()`. `mutation_length` and `clock_length` require partition access. Other fields come from `NodeTimetree`, `EdgeTimetree`, `ClockModel`, `NodeConfidenceInterval`.

## Dependencies

- `io-implement-augur-node-data-json-types.md` (serde types)
- `dates-preserve-raw-string-and-input-type.md` (raw_date, date, date_inferred fields)

## Reference

- Format specification: [../reports/augur-node-data-json.md](../reports/augur-node-data-json.md) section "`augur refine` node data JSON"
- Augur transformation code: [`refine.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/refine.py)
