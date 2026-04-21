# AlphabetConfig validation misses unknown-inside-ambiguous check

`AlphabetConfig::validate()` at [packages/treetime/src/alphabet/alphabet_config.rs#L137-L145](../../packages/treetime/src/alphabet/alphabet_config.rs#L137-L145) contains two consecutive blocks that both test `ambiguous_keys.contains(*gap)`. The second block (line 142) was intended to check `ambiguous_keys.contains(*unknown)` but instead duplicates the gap check. The error message on line 144 reads "Ambiguous set contains 'unknown' character" but the condition tests for `gap`.

An alphabet configuration where the `unknown` character appears as an ambiguous key passes validation and produces a malformed `Alphabet` where `set_to_char` maps a single character to two distinct roles (unknown and ambiguous).

## Affected code

- Duplicate check: [packages/treetime/src/alphabet/alphabet_config.rs#L142](../../packages/treetime/src/alphabet/alphabet_config.rs#L142) -- `ambiguous_keys.contains(*gap)` should be `ambiguous_keys.contains(*unknown)`
- Existing gap check (correct): [packages/treetime/src/alphabet/alphabet_config.rs#L137](../../packages/treetime/src/alphabet/alphabet_config.rs#L137)

## Fix

Change `ambiguous_keys.contains(*gap)` on line 142 to `ambiguous_keys.contains(*unknown)`.

Add a unit test `test_alphabet_config_validate_ambiguous_contains_unknown` mirroring the existing `test_alphabet_config_validate_ambiguous_contains_gap` test.
