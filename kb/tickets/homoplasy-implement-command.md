# Implement homoplasy command

`run_homoplasy()` at [packages/treetime/src/commands/homoplasy/run.rs#L5](../../packages/treetime/src/commands/homoplasy/run.rs#L5) unconditionally panics with `unimplemented!()`. The CLI accepts `treetime homoplasy` and parses all flags via `TreetimeHomoplasyArgs` at [packages/treetime/src/commands/homoplasy/args.rs](../../packages/treetime/src/commands/homoplasy/args.rs), but the command produces a panic rather than output.

## Scope

v0 homoplasy command computes and reports homoplasies (recurrent mutations) on a phylogenetic tree. The v0 implementation includes:

- Identification of homoplastic mutations by comparing ancestral reconstruction to observed sequences
- Detailed reporting of per-site homoplasy counts
- DRM (drug resistance mutation) annotation
- Rescaling of branch lengths by homoplasy density

None of these features are implemented in v1.

## Affected code

- Panic site: [packages/treetime/src/commands/homoplasy/run.rs#L5](../../packages/treetime/src/commands/homoplasy/run.rs#L5)
- CLI args (parsed but unused): [packages/treetime/src/commands/homoplasy/args.rs](../../packages/treetime/src/commands/homoplasy/args.rs)
- CLI dispatch: [packages/treetime-cli/src/bin/treetime.rs](../../packages/treetime-cli/src/bin/treetime.rs)

## Fix

Replace `unimplemented!()` with a descriptive error message: "The homoplasy command is not yet implemented in v1." Implement the command or remove it from the CLI.

## Related issues

- Source: [H-homoplasy-command-unimplemented.md](../issues/H-homoplasy-command-unimplemented.md) -- delete after full resolution
- [unimplemented.md](../algo/unimplemented.md) -- homoplasy algorithms listed as unimplemented
- [README.md](../features/README.md) -- homoplasy features tracked as 0/9 ported
