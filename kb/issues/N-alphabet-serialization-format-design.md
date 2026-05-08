# Alphabet serialization format design

## Summary

`Alphabet` serialization needs a user-facing format design before release. The current implementation routes serde through `AlphabetConfig`, which serializes raw byte values (`u8`, `Vec<u8>`) instead of human-readable characters. No command currently outputs alphabet information, but alphabet output will be needed when features like automatic alphabet deduction from data are introduced.

## Current state

`Alphabet` stores an `AlphabetConfig` and uses `#[serde(try_from = "AlphabetConfig")]` for deserialization and delegates `Serialize` to the stored config. `AlphabetConfig` fields are raw bytes:

- `canonical: Vec<u8>` - serializes as `[65, 67, 71, 84]` instead of `["A", "C", "G", "T"]`
- `ambiguous: IndexMap<u8, Vec<u8>>` - keys and values are byte codes
- `unknown: u8`, `gap: u8` - single byte codes

No production code serializes `Alphabet` today. The `Serialize` impl exists from the original `derive(Serialize, Deserialize)` and is only exercised by roundtrip tests.

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

- Keep `AlphabetConfig` as `u8` internally, add serde attributes for human-readable serialization
- Change `AlphabetConfig` fields to `AsciiChar` / `Vec<AsciiChar>` throughout
- Separate internal config from a public-facing schema type

## Locations

- `packages/treetime/src/alphabet/alphabet.rs` - `struct Alphabet`, `Serialize` impl, `TryFrom<AlphabetConfig>`
- `packages/treetime/src/alphabet/alphabet_config.rs` - `struct AlphabetConfig`
- `packages/treetime/src/commands/*/args.rs` - `--alphabet` CLI flag (4 commands)
- `packages/treetime-io/src/fasta.rs` - `FastaReader` uses `AlphabetLike` trait

## Impact

Low priority. No production serialization exists. Blocking factor for any feature that reports alphabet information back to the user.
