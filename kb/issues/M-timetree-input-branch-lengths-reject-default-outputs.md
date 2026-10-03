# Timetree with input branch lengths fails with the default outputs

`treetime timetree --branch-length-mode=input --output-all=<dir>` always fails before inference:

```
Error:
   0: Reconstructed sequence output requires ancestral reconstruction; incompatible with --branch-length-mode=input
```

The default selection of `--output-all` includes the reconstructed nucleotide FASTA [packages/app-output/src/output_plan.rs#L240-L256](../../packages/app-output/src/output_plan.rs#L240-L256), and `fn validate_params` rejects any sequence output in input mode [packages/treetime/src/timetree/pipeline.rs#L229-L240](../../packages/treetime/src/timetree/pipeline.rs#L229-L240). A user therefore has to know to pass an `--output-selection` without the FASTA. The failed run leaves an empty `timetree.reconstructed-nuc.fasta` in the output directory.

## Open question

What should input mode do with the default outputs?

- Leave the reconstructed FASTA out of the default selection in input mode, and keep the error for an explicit request
- Skip the FASTA with a warning, as `--output-all` already does for confidence intervals that were not computed

## Smoke coverage

The `branch-input-all-outputs` row of `dev/smoke.toml` declares this failure. The `branch-input` rows select outputs without the FASTA, so they reach the time inference, where [H-timetree-input-branch-lengths-abort-on-point-division.md](H-timetree-input-branch-lengths-abort-on-point-division.md) stops them.
