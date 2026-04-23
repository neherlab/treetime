# Ancestral joint method default panics

`MethodAncestral::Joint` is the `#[default]` variant of the CLI `--method-anc` enum, so running `treetime ancestral` without an explicit `--method-anc` flag dispatches to the `Joint` branch at [packages/treetime/src/commands/ancestral/run.rs#L198](../../packages/treetime/src/commands/ancestral/run.rs#L198), which unconditionally panics with `unimplemented!()`.

Joint reconstruction was intentionally removed from v1 (documented in [port-intentional-changes/ancestral-joint-reconstruction-removed.md](../port-intentional-changes/ancestral-joint-reconstruction-removed.md)), but the `#[default]` annotation on the `Joint` variant was not updated to reflect this removal. The result is that the most common invocation of the `ancestral` command crashes by construction.

## Affected code

- Enum definition: [packages/treetime/src/commands/ancestral/args.rs#L11-L16](../../packages/treetime/src/commands/ancestral/args.rs#L11-L16) -- `#[default]` on `Joint`
- Panic site: [packages/treetime/src/commands/ancestral/run.rs#L199](../../packages/treetime/src/commands/ancestral/run.rs#L199) -- `unimplemented!()`

## Reproduction

```bash
treetime ancestral --tree=data/flu/h3n2/20/tree.nwk data/flu/h3n2/20/aln.fasta.xz
```

Panics with `not implemented` at `run.rs:199`.

## Fix

Move `#[default]` from `Joint` to `Marginal`. Optionally remove the `Joint` variant entirely, or replace the `unimplemented!()` with a descriptive error message directing users to `--method-anc=marginal` or `--method-anc=parsimony`.

## Related

- [port-intentional-changes/ancestral-joint-reconstruction-removed.md](../port-intentional-changes/ancestral-joint-reconstruction-removed.md): documents the rationale for removing joint reconstruction
- [port-algo-inventory/ancestral.md](../port-algo-inventory/ancestral.md): algorithm inventory
