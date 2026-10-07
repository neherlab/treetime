# Alphabet serialization format design

## Summary

Reading or writing an alphabet needs a user-facing format design before release. `Alphabet` and `AlphabetConfig` have no serialization: no command reads a custom alphabet or outputs alphabet information. Alphabet output will be needed when features like automatic alphabet deduction from data are introduced, and alphabet input when users can define custom alphabets.

## Current state

`Alphabet::with_config()` builds an alphabet from an `AlphabetConfig`, and `Alphabet::new()` builds the three named alphabets from fixed configs. `AlphabetConfig` fields are raw bytes:

- `canonical: Vec<u8>` - byte codes such as `[65, 67, 71, 84]` instead of `["A", "C", "G", "T"]`
- `ambiguous: IndexMap<u8, Vec<u8>>` - keys and values are byte codes
- `unknown: u8`, `gap: u8` - single byte codes

A direct derive on these fields would read and write byte codes (or, for byte vectors, base64 strings), which a user cannot read or write by hand.

## Design scope

### Output (not yet implemented)

Commands do not emit alphabet information. Anticipated needs:

- Reporting the deduced alphabet when auto-detected from input data
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
