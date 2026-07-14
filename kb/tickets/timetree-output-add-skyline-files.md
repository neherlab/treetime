# Add skyline output files

v0 `--coalescent=skyline` produces `skyline.tsv` (effective population size over time with confidence bounds) and `skyline.pdf`. v1 computes the `SkylineResult` containing `time_grid`, `log_tc_values`, `log_likelihood` in memory but never serializes it to disk.

## Location

`SkylineResult` computed at [`run.rs#L285-L291`](../../packages/treetime/src/commands/timetree/run.rs#L285-L291) but not written.

## Details

- v0 `skyline.tsv` columns: `#date`, `N_e`, `lower`, `upper` (confidence bounds)
- v0 defaults to 20 skyline grid points, v1 defaults to 10
  ([`args.rs#L175`](../../packages/treetime/src/commands/timetree/args.rs#L175))

## Related issues

- Source: [N-timetree-missing-skyline-output.md](../issues/N-timetree-missing-skyline-output.md) -- delete after full resolution
- [kb/issues/N-timetree-missing-output-files.md](../issues/N-timetree-missing-output-files.md)
  umbrella issue for all missing outputs
- [timetree.json missing coalescent parameters](timetree-output-json-missing-coalescent-parameters.md)
  skyline results also not serialized to JSON
- [Skyline coalescent uses Nelder-Mead instead of SLSQP](../issues/M-timetree-skyline-nelder-mead-optimizer.md)
  different optimizer produces different skyline values
