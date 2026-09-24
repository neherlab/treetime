# Command config schemas carry the description of a flattened argument group

The root `description` of each `packages/schemas/input-config-<command>.schema.json` is the doc comment of the first flattened argument group, not a description of the command. For `timetree`, `clock`, `ancestral`, `optimize` and `prune` it is the text of `struct AlignmentArgs` ("Sequence alignment input shared by all commands that read sequences..."); for `mugration` it is the text of `struct MetadataIdArgs`. The `Treetime<Command>ArgsRaw` structs have no doc comment of their own, and the generated root description is that of their first flattened field type.

The same text appears as the description of the `<Command>Config` components of `packages/app-contracts/openapi.json` and of the generated TypeScript types.

## Impact

Editors that show the schema description for a config file, and a settings form that shows it as the command description, show text about alignments or metadata columns.

## Fix direction

Give each `Treetime<Command>ArgsRaw` a schema description of its own, for example with `#[schemars(description = "...")]` taken from the command's `about` text, and regenerate the schemas and the TypeScript client.
