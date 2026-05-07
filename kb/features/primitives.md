# Sequence Primitives

- [x] AsciiChar (validated u8, type-safe character)
- [x] Seq (Vec\<AsciiChar\> with string-like operations, concatenation, repetition)
- [x] BitSet128/StateSet (u128 bitfield for character sets, set operations)
- [x] Sub (substitution mutation: pos + ref + qry, invertible, parseable)
- [x] InDel (insertion/deletion: range + seq + direction, invertible)
- [x] Composition (character frequency tracking with mutation updates)
- [x] find_char_ranges (contiguous range detection for gaps, unknowns)
- [x] compute_divs (root-to-tip divergence calculation)
