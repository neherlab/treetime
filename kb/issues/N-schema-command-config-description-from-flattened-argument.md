# Command config schemas carry the description of a flattened argument group

The root `description` of each `packages/schemas/input-config-<command>.schema.json` is the doc comment of the first flattened argument group, not a description of the command:

- `timetree`, `clock`, `ancestral`, `optimize`, `prune` and `homoplasy`: the text of `struct AlignmentArgs` ("Sequence alignment input shared by all commands that accept an alignment...") [packages/app-commands/src/commands/shared/alignment.rs#L21-L32](../../packages/app-commands/src/commands/shared/alignment.rs#L21-L32). The `clock` schema thus describes an alignment input that `clock` does not read ([M-clock-alignment-ignored.md](M-clock-alignment-ignored.md))
- `mugration`: the text of `struct TreeDialectArgs` ("Newick dialect of the tree inputs, shared by every command that reads a tree.") [packages/app-commands/src/commands/shared/tree_input.rs#L40-L48](../../packages/app-commands/src/commands/shared/tree_input.rs#L40-L48)

The `Treetime<Command>ArgsRaw` structs have no doc comment of their own, so the generated root description is that of their first flattened field type, for example `alignment` in `struct TreetimeTimetreeArgsRaw` [packages/app-commands/src/commands/timetree/args.rs#L113-L123](../../packages/app-commands/src/commands/timetree/args.rs#L113-L123).

The same text appears as the description of the `<Command>Config` components of `packages/app-contracts/openapi.json` and in the generated TypeScript types and zod schemas (`packages/app-contracts/src/generated/types.gen.ts`, `packages/app-contracts/src/generated/zod.gen.ts`).

## Impact

Editors that show the schema description for a config file, and a settings form that shows it as the command description, show text about the alignment input or the Newick dialect.

## Fix direction

Give each `Treetime<Command>ArgsRaw` a schema description of its own, for example with `#[schemars(description = "...")]` taken from the command's `about` text, and regenerate the schemas and the TypeScript client.
