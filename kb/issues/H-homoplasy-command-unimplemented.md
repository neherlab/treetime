# Homoplasy command is unimplemented

`fn run_homoplasy()` [`packages/treetime/src/commands/homoplasy/run.rs#L6-L8`](../../packages/treetime/src/commands/homoplasy/run.rs#L6-L8) returns an explicit not-implemented error. The CLI accepts `treetime homoplasy` and parses all flags via `struct TreetimeHomoplasyArgs` [`packages/treetime/src/commands/homoplasy/args.rs#L12-L28`](../../packages/treetime/src/commands/homoplasy/args.rs#L12-L28), but the command cannot produce homoplasy results.

## Scope

v0 homoplasy command computes and reports homoplasies (recurrent mutations) on a phylogenetic tree. The v0 implementation includes:

- Identification of homoplastic mutations by comparing ancestral reconstruction to observed sequences
- Detailed reporting controlled by a flag
- DRM (drug resistance mutation) annotation
- Uniform multiplication of branch lengths by the numeric `--rescale` factor

None of these features are implemented in v1.

## Affected code

- Error site: [`packages/treetime/src/commands/homoplasy/run.rs#L6-L8`](../../packages/treetime/src/commands/homoplasy/run.rs#L6-L8)
- CLI args (parsed but unused): [`packages/treetime/src/commands/homoplasy/args.rs#L12-L28`](../../packages/treetime/src/commands/homoplasy/args.rs#L12-L28)
- CLI dispatch: [`packages/app-cli/src/bin/treetime.rs#L104-L105`](../../packages/app-cli/src/bin/treetime.rs#L104-L105)

## Decisions required

The current v1 argument types do not match v0: v1 parses `--rescale` as a boolean and `--detailed` as an optional string, while v0 uses a numeric rescaling factor and a detail flag. [`packages/treetime/src/commands/homoplasy/args.rs#L20-L26`](../../packages/treetime/src/commands/homoplasy/args.rs#L20-L26) [`packages/legacy/treetime/treetime/argument_parser.py#L372-L376`](../../packages/legacy/treetime/treetime/argument_parser.py#L372-L376)

Before creating an implementation ticket, specify exact CLI parity, output schemas, numerical and statistical oracles, VCF behavior, DRM validation, and the complete acceptance contract. The existing explicit error is the correct failure mode while the command is unavailable.

## Related

- [algo/unimplemented.md](../algo/unimplemented.md): homoplasy algorithms listed as unimplemented
- [features/\README.md](../features/README.md): homoplasy features tracked as 0/9 ported
