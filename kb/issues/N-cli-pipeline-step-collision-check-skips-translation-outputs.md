# Pipeline step collision check skips translation outputs

`treetime pipeline` checks the plan before any step runs and rejects two steps that write the same file (`fn reject_output_collisions()`, [packages/app-cli/src/cli/pipeline/safety.rs#L56-L71](../../packages/app-cli/src/cli/pipeline/safety.rs#L56-L71)). The check reads the outputs of each step from `fn output_paths()` ([packages/app-cli/src/cli/pipeline/safety.rs#L104-L113](../../packages/app-cli/src/cli/pipeline/safety.rs#L104-L113)), which drops every path whose file name holds the `{cds}` placeholder: the per-CDS amino-acid FASTA outputs of `--output-translations`. Two steps whose translation outputs expand to the same file therefore both run, and the later step overwrites the file of the earlier one.

The single-command check expands these templates: `RunPaths::new()` keeps the translation outputs apart as `cds_outputs`, and `RunPaths::check_with_cds_names()` expands them with the names of `--cds-selection` and compares them with the inputs of the run ([packages/app-commands/src/commands/shared/run_paths.rs#L87-L136](../../packages/app-commands/src/commands/shared/run_paths.rs#L87-L136)). The pipeline check has no such expansion.

## Example

Two `ancestral` steps on `data/rsv/a/20` with `translations: ".../translations/{cds}.fasta.xz"` and `cds_selection: ["F"]`, both with `output_translations: "aa/{cds}.fasta"`, pass the plan check: the run exits 0, and both steps write `aa/F.fasta` in the folder of the config file. The same two steps with `output_sequences: "aa/F.fasta"` stop before the first step:

```text
Error:
   0: steps `a` and `b` both write `<config folder>/aa/F.fasta`; each output path must be written by only one step
```

A translation output of one step that expands to the plain output of another (step `a` with `output_translations: "aa/{cds}.fasta"` and step `b` with `output_sequences: "aa/F.fasta"`) also passes, by the same filter.

The same `fn output_paths()` feeds `fn reject_missing_upstream()` ([#L73-L102](../../packages/app-cli/src/cli/pipeline/safety.rs#L73-L102)), so a selected step that reads a translation file of an unselected step is not checked for an absent upstream output either.

## Fix direction

Expand the translation templates of each step with its `cds_selection` names, as `RunPaths::check_with_cds_names()` does, before the collision and upstream checks. A step without `cds_selection` takes its CDS names from the annotation at run time; for such a step, compare the templates themselves (two equal templates collide, and a template collides with a plain path that it can expand to), or read the CDS names of the annotation during the plan check. Add a case to the `safety.rs` tests for each collision kind above.
