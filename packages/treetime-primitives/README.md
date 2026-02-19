# treetime-primitives

Primitive types for sequence representation in treetime. Provides the foundational character, sequence, and bitset types used throughout the treetime crate ecosystem.

## Types

### `AsciiChar`

Single ASCII character wrapper around `u8`. Used as the element type for sequences and bitset membership. Supports conversions to/from `u8`, `u16`, `u32`, `u64`, `usize`, and `char`.

```rust
use treetime_primitives::AsciiChar;

let c = AsciiChar::try_new(b'A')?;  // fallible
let c = AsciiChar::from_byte_unchecked(b'A');  // unchecked for known-valid literals
let n: u8 = c.into();
```

### `Seq`

Ordered collection of `AsciiChar` values representing a genetic sequence. Backed by `Vec<AsciiChar>` with a `Vec`-like API (push, pop, indexing, slicing, iteration). Supports zero-copy `as_str()` conversion, `serde` serialization, and `std::io::Read`/`Write`.

```rust
use treetime_primitives::Seq;

let seq = Seq::try_from_str("ACGT")?;
assert_eq!(seq.len(), 4);
assert_eq!(seq.as_str(), "ACGT");
```

### `BitSet128`

Compact set of up to 128 elements stored as a `u128` bitmask. Supports set operations (union, intersection, difference, symmetric difference) via methods and operators (`|`, `&`, `-`, `^`). Collects from iterators of `AsciiChar`, `u8`, or `char`.

```rust
use treetime_primitives::{BitSet128, bitset128};

let set = bitset128!['A', 'C', 'G', 'T'];
assert_eq!(set.len(), 4);
assert!(set.contains(b'A'));

let other = bitset128!['A', 'G'];
let union = set | other;
let intersection = set & other;
```

`BitSet128Status` classifies a bitset as `Empty`, `Unambiguous(AsciiChar)`, or `Ambiguous(BitSet128)`.

### `StateSet` / `stateset!`

Type alias for `BitSet128`, used to represent possible character states at a sequence position (e.g., ambiguous nucleotides). The `stateset!` macro mirrors `bitset128!`.

### `AlphabetLike`

Trait for types that can validate characters and enumerate their character set. Decouples I/O operations from the full `Alphabet` type in the core crate.

```rust
use treetime_primitives::{AlphabetLike, AsciiChar};

fn validate<A: AlphabetLike>(alphabet: &A, c: AsciiChar) -> bool {
    alphabet.contains(c)
}
```

## Dependents

- `treetime` - core library
- `treetime-io` - file I/O
