# Alphabet serialization format design

## Summary

Reading or writing an alphabet needs a user-facing format design before release. `Alphabet` and `AlphabetConfig` have no serialization: no command reads a custom alphabet or outputs alphabet information. Alphabet output is needed to report the alphabet that a run used, and alphabet input when users can define custom alphabets.

> [!IMPORTANT]
> **Decision required.** The user-facing alphabet format and its internal representation are not chosen. Options for the representation are listed under "Internal representation"; the format needs the choices listed under "Format considerations". No command needs the format yet, so no option is forced by current code.

## Current state

`Alphabet::with_config()` [packages/treetime/src/alphabet/alphabet.rs#L102](../../packages/treetime/src/alphabet/alphabet.rs#L102) builds an alphabet from an `AlphabetConfig`, and `Alphabet::new()` [packages/treetime/src/alphabet/alphabet.rs#L47](../../packages/treetime/src/alphabet/alphabet.rs#L47) builds the three named alphabets from fixed configs. `AlphabetConfig` [packages/treetime/src/alphabet/alphabet_config.rs#L13-L18](../../packages/treetime/src/alphabet/alphabet_config.rs#L13-L18) fields are raw bytes:

- `canonical: Vec<u8>` - byte codes such as `[65, 67, 71, 84]` instead of `["A", "C", "G", "T"]`
- `ambiguous: IndexMap<u8, Vec<u8>>` - keys and values are byte codes
- `unknown: u8`, `gap: u8` - single byte codes

A direct derive on these fields would read and write byte codes (or, for byte vectors, base64 strings), which a user cannot read or write by hand.

The sequence commands do not deduce the alphabet from the data. When `--alphabet` is omitted, every caller takes `AlphabetName::default()`, which is `nuc`, through `self.alphabet.unwrap_or_default()` [packages/app-commands/src/commands/shared/sequence_inputs.rs#L105](../../packages/app-commands/src/commands/shared/sequence_inputs.rs#L105).

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

## Proposed implementation plan

The internal representation and the format stay open (Decision required above). Common steps first, then one short plan per representation option.

Common steps:

1. Define the user-facing schema: `canonical` as a list of one-character strings, `ambiguous` as a map from one character to a list of characters, `gap` and `unknown` as one character each. Derive `Serialize`, `Deserialize`, and `JsonSchema` on the schema type, so `just gen` writes its JSON schema into `packages/schemas`
2. Validate at the boundary: every character is one ASCII byte, `gap` and `unknown` are not canonical, every ambiguity code maps to canonical characters only. The checks of `fn Alphabet::with_config()` (`packages/treetime/src/alphabet/alphabet.rs`) and `fn AlphabetConfig.create_profile_map()` (`packages/treetime/src/alphabet/alphabet_config.rs`) stay the single owner; the reader calls them and returns their errors with the file path in context
3. Output: add the schema value of the alphabet in use to the structured output that the first consumer needs (run report or node-data JSON), at the place that writes that output
4. Input: when a command needs custom alphabets, extend `--alphabet` in `packages/app-commands/src/commands/shared/alphabet.rs` to accept a file path besides the named alphabets (`nuc`, `aa`, `aa-no-stop`)

Option plans:

- **Keep `u8`, add adapters**: field-level `#[serde(with = ...)]` helpers on `struct AlphabetConfig` that read and write characters. Smallest change; the byte fields stay
- **Change fields to `AsciiChar`**: change `canonical`, `ambiguous`, `unknown`, `gap` to `AsciiChar` types and update `fn Alphabet::with_config()` and `fn AlphabetConfig.create_profile_map()`. Serialization of `AsciiChar` as a one-character string then serves the schema directly
- **Separate schema type**: a public schema struct with `String`/`char` fields and a parsing conversion into `AlphabetConfig`; `AlphabetConfig` stays internal

Tests (in `packages/treetime/src/alphabet/__tests__/`):

- Round-trip property test: schema value of each named alphabet, written and read again, builds an `Alphabet` equal to `Alphabet::new(name)` (oracle: the fixed configs in `fn Alphabet::new()`)
- Unit tests with `indoc!` JSON for the `nuc` alphabet (oracle: the IUPAC nucleotide code table), and one `assert_err!` case per validation rule of step 2

Validation: `just check`, `just test-rs`, `just generated-check` after `just gen`.
