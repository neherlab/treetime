# Smoke reports config cases as changed because outputs hold the snapshot path

## Problem

A config case of `dev/smoke` runs `treetime <command> --config=<file>` in a sandbox directory inside the snapshot directory ([dev/smoke.toml](../../dev/smoke.toml), section `[config]`). The config reader resolves the relative input paths of the file to absolute paths, and the augur node data output records them in its `alignment` and `input_tree` fields. These paths contain the snapshot directory, which differs between the baseline and the current snapshot.

`fn normalized_stderr` in [dev/smoke](../../dev/smoke) replaces the snapshot directory with `<snapshot>`, but `fn normalized_bytes` removes only the fields listed in `VOLATILE_FIELDS`. The output files therefore differ in every comparison, and the case is reported as `changed` even when the command output is the same.

## Evidence

`./dev/smoke --only optimize` on a branch with identical `optimize` results reported one `changed` case, `config/flu/h3n2/20/optimize`, with this difference in `optimize.augur-node-data.json`:

```
alignment:  ".../snapshots/<baseline>/cases/config/flu/h3n2/20/optimize/sandbox/data/flu/h3n2/20/aln.fasta.xz" -> ".../snapshots/<current>/cases/..."
input_tree: ".../snapshots/<baseline>/cases/config/flu/h3n2/20/optimize/sandbox/data/flu/h3n2/20/tree.nwk"     -> ".../snapshots/<current>/cases/..."
```

## Fix direction

Replace the snapshot directory with `<snapshot>` in the compared output files, as `fn normalized_stderr` does for stderr.

## Validation

A smoke comparison of two snapshots with the same CLI outputs reports every config case as identical.
