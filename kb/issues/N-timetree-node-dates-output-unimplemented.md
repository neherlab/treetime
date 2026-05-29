# write_node_dates() is a todo!() stub

`write_node_dates()` is unimplemented. No `dates.tsv` output file is produced. v0 writes `dates.tsv` with columns including node name, date estimate, and (with `--confidence --covariation`) lower/upper bounds.

The same node date data (node name, `numdate`, resolved `date`, and `num_date_confidence` lower/upper bounds) is now emitted in `timetree.augur-node-data.json`. The tabular `dates.tsv` file and the `write_node_dates()` stub remain unimplemented.

## Related issues

- [Missing output files compared to v0](N-timetree-missing-output-files.md)
  umbrella issue for all missing outputs
