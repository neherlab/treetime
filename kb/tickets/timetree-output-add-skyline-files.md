# Add skyline output files

Let $T_c$ denote the coalescent population-size time scale. v0 `--coalescent=skyline` produces `skyline.tsv` (effective population size over time with confidence bounds) and `skyline.pdf`. v1 computes a `SkylineResult`, retains only its $T_c$ distribution for the final inference pass, and never serializes the result to disk.

## Location

`SkylineResult` is computed at [`packages/treetime/src/timetree/pipeline.rs#L343-L351`](../../packages/treetime/src/timetree/pipeline.rs#L343-L351) but not returned to the command output or written.

## Details

- v0 `skyline.tsv` columns: `#date`, `N_e`, `lower`, `upper` (confidence bounds)
- v0 defaults to 20 skyline grid points, v1 defaults to 10
  ([`packages/treetime/src/commands/timetree/args.rs#L142-L149`](../../packages/treetime/src/commands/timetree/args.rs#L142-L149))

## Related issues

- Source: [N-timetree-missing-skyline-output.md](../issues/N-timetree-missing-skyline-output.md) -- delete after full resolution
- [kb/issues/N-timetree-missing-output-files.md](../issues/N-timetree-missing-output-files.md)
  umbrella issue for all missing outputs
- [timetree.json missing coalescent parameters](timetree-output-json-missing-coalescent-parameters.md)
  skyline results also not serialized to JSON
- [Skyline coalescent uses Nelder-Mead instead of SLSQP](../issues/M-timetree-skyline-nelder-mead-optimizer.md)
  different optimizer produces different skyline values
