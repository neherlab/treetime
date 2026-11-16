# CLI run warnings and diagnostics reach stderr only

The CLI writes warnings only to the log on stderr (`fn write_log()` in [packages/app-cli/src/cli/progress.rs](../../packages/app-cli/src/cli/progress.rs)). A pipeline that keeps only the output files loses them, and no output file records which inputs, settings and version produced a result. The 31 `progress_warn!` sites outside tests are untyped text.

The app already keeps a typed run record: `struct RunRecord` with version, command, inputs with SHA-256, `config_hash`, warnings and error ([packages/app-backend/src/runs/record.rs](../../packages/app-backend/src/runs/record.rs)), filled by `struct WarningCollector` ([packages/app-backend/src/runs/warnings.rs#L4](../../packages/app-backend/src/runs/warnings.rs#L4)).

Codon decoding and the annotation reader need a place for aggregate diagnostics per CDS: diffuse codons, sites where decoding rules differ, decoded internal stop codons and stop mass, overlap conflicts, partly gapped internal codons, strands read as forward, gene fallbacks, phase mismatches, the coordinate check, and consumer checks (every child mutation starts from its parent's state).

## Fix direction

One shared typed run record for the CLI and the app ([kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md), axis 13): move the collector out of `app-backend` into a shared crate, give the CLI a record file with the fields of `RunRecord`, add a `diagnostics` struct with one optional named section per subsystem, and turn the untyped warning sites into `RunWarningKind` variants case by case. Per-site diagnostic lists stay outside the first version; they need their own format.

> [!IMPORTANT]
> **Decision required.** The name `report.json` collides with the homoplasy report (`--output-homoplasy-report`, `.report.txt`) and the option `--report-ambiguous`. Options: `report.json`, keeping the qualified name "homoplasy report"; or another name, such as `run.json`. Also open: whether the record carries expected substitution counts in addition to counts from point estimates.

## Validation

- A CLI run that raises a warning writes the record with that warning as a typed entry, and the app shows the same entry for the same run
- The record schema is generated with the other schemas and checked by `just generated-check`
