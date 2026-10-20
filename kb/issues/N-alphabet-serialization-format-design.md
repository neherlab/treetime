# Alphabet serialization format design

## Summary

Reading or writing an alphabet needs a user-facing format design before release. `Alphabet` and `AlphabetConfig` have no serialization: no command reads a custom alphabet or outputs alphabet information. Alphabet output is needed to report the alphabet that a run used, and alphabet input when users can define custom alphabets.

> [!IMPORTANT]
> **Decision required.** The user-facing alphabet format and its internal representation are not chosen. Options for the representation are listed under "Internal representation"; the format needs the choices listed under "Format considerations". No command needs the format yet, so no option is forced by current code.

## Current state

`Alphabet::with_config()` [packages/treetime/src/alphabet/alphabet.rs#L89](../../packages/treetime/src/alphabet/alphabet.rs#L89) builds an alphabet from an `AlphabetConfig`, and `Alphabet::new()` [packages/treetime/src/alphabet/alphabet.rs#L43](../../packages/treetime/src/alphabet/alphabet.rs#L43) builds the three named alphabets from fixed configs. `AlphabetConfig` [packages/treetime/src/alphabet/alphabet_config.rs#L13-L18](../../packages/treetime/src/alphabet/alphabet_config.rs#L13-L18) fields are raw bytes:

- `canonical: Vec<u8>` - byte codes such as `[65, 67, 71, 84]` instead of `["A", "C", "G", "T"]`
- `ambiguous: IndexMap<u8, Vec<u8>>` - keys and values are byte codes
- `unknown: u8`, `gap: u8` - single byte codes

A direct derive on these fields would read and write byte codes (or, for byte vectors, base64 strings), which a user cannot read or write by hand.

The sequence commands do not deduce the alphabet from the data. When `--alphabet` is omitted, every caller takes `AlphabetName::default()`, which is `nuc`, through `alphabet_name().unwrap_or_default()`, for example [packages/app-commands/src/commands/shared/sequence_inputs.rs#L26](../../packages/app-commands/src/commands/shared/sequence_inputs.rs#L26).

## Design scope

### Output (not yet implemented)

Commands do not emit alphabet information. Anticipated needs:

- Reporting the alphabet that a run used, in particular when the user did not select one
- Including alphabet metadata in structured output files (JSON, Nexus annotations)
- Diagnostic output for debugging character-set mismatches

Format considerations:

- Human readability: characters (`"A"`) not byte codes (`65`)
- Structured output: canonical set, ambiguity mappings, gap/unknown markers
- Consistency with other output formats in the codebase (GTR JSON, auspice JSON)

### Input

`AlphabetConfig` is the natural input format for custom alphabets (canonical characters, ambiguity mappings, gap, unknown). Currently alphabets are selected via `--alphabet` CLI flag from a fixed enum (`nuc`, `aa`, `aa-no-stop`). Custom alphabet configs from files would use `AlphabetConfig` deserialization, which needs a human-writable format (characters, not byte codes).

### Internal representation

`AlphabetConfig` stores `u8` because it predates `AsciiChar`. Options:

- Keep `AlphabetConfig` as `u8` internally, add adapters for human-readable serialization
- Change `AlphabetConfig` fields to `AsciiChar` / `Vec<AsciiChar>` throughout
- Separate internal config from a public-facing schema type

## Locations

- `packages/treetime/src/alphabet/alphabet.rs` - `struct Alphabet`, `fn Alphabet::with_config()`
- `packages/treetime/src/alphabet/alphabet_config.rs` - `struct AlphabetConfig`
- `packages/app-commands/src/commands/shared/alphabet.rs` - `--alphabet` CLI flag shared by the sequence commands
- `packages/treetime-io/src/fasta.rs` - `fasta_read()` uses `AlphabetLike` trait

## Impact

Low priority. No alphabet is read from or written to a file. Blocking factor for any feature that reports alphabet information back to the user.
