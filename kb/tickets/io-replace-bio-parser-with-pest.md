# Replace bio Newick parser with pest

Replace `bio::io::newick::read` with a custom pest-based parser for base Newick (no comment parsing yet). Parse branch lengths as f64 directly. Remove the bio crate dependency.

## Scope

- Write a pest grammar for base Newick (~20 lines): Tree, Subtree, Leaf, Internal, BranchSet, Branch, Name (quoted/unquoted), Length
- Implement parser producing `Graph<N, E, D>` via existing `NodeFromNwk`/`EdgeFromNwk` traits
- Parse branch lengths as f64 (bio parses f32 at `bio-2.3.0/newick.rs:106`, widened to f64 too late at `nwk.rs:80`)
- Strip `[...]` content as plain comments (same behavior as bio, but without erroring)
- Pass empty `BTreeMap` to `NodeFromNwk::from_nwk` (same as current `nwk.rs:73` -- comment parsing is ticket #3)
- Add `pest` and `pest_derive` as direct workspace dependencies
- Remove `bio` crate dependency (verify no other uses first)
- Preserve existing `nwk_read`, `nwk_read_file`, `nwk_read_str` API signatures

## Grammar reference

Base Newick grammar from [kb/reports/newick-annotation-dialects.md](../reports/newick-annotation-dialects.md) lines 23-32. Bio's existing pest grammar at `bio-2.3.0/src/io/newick.pest` (19 lines) for comparison -- its `safe` rule at line 11 excludes `[` and `]`, which is the root cause of the parse error.

## Tests

### Happy paths -- tree structure

- Binary tree: `(A:0.1,B:0.2);`
- Named internal: `(A:0.1,B:0.2)root;`
- No branch lengths: `(A,B);`
- Mixed lengths (some present, some absent): `(A:0.1,B,(C:0.3,D));`
- Deep nesting: `((((A,B),C),D),E);`
- Multifurcating: `(A,B,C,D);`
- Root with branch length: `(A:0.1,B:0.2):0.0;`
- Quoted names: `('sample name':0.1,'other name':0.2);`
- Underscore-to-space mapping in unquoted labels: `Sample_1` -> `Sample 1`

### Happy paths -- branch lengths

- Standard decimal: `0.1`, `0.005`, `123.456`
- Scientific notation: `1.5e-3`, `2.0E-5`, `1e10`
- Very small: `1e-15` (precision preserved)
- Very large: `1e5`
- Zero: `0.0`, `0`
- f64 precision: `0.005` reads as `0.005000000000000000` not `0.004999999888241291` (f32 artifact)
- Many decimal places: `0.123456789012345` (15 significant digits preserved)

### Happy paths -- comments

- Comment after name: `(A[comment]:0.1,B:0.2);` -> comment stripped, tree parsed
- Comment after length: `(A:0.1[comment],B:0.2);` -> stripped
- BEAST-style comment (stripped, not parsed -- ticket #3): `(A[&mutations="G42T"]:0.1,B:0.2);`
- Multiple comments on one node: `(A[c1][c2]:0.1,B:0.2);`
- Comment on internal node: `(A,B)[comment];`

### Edge cases

- Single leaf: `A;`
- Single leaf with length: `A:0.1;`
- Empty name (unnamed internal node): `(:0.1,:0.2);`
- Name parseable as number: `(42:0.1,B:0.2);` -- current bio behavior: discards as weight. Preserve this behavior
- Whitespace everywhere: `( A : 0.1 , B : 0.2 ) root ;`
- No whitespace at all: `(A:0.1,B:0.2)root;`
- Trailing newlines and whitespace after semicolon
- Tab characters in whitespace positions
- Quoted name containing single quotes: `('it''s a name':0.1)` (Newick escapes `'` as `''`)
- Quoted name containing parentheses: `('(A)':0.1,B:0.2);`
- Negative branch length: `(A:-0.001,B:0.2);` (some tools emit these)
- Branch length with leading dot: `(A:.1,B:.2);`
- Quoted empty name: `('':0.1,B:0.2);`

### Pathological -- error handling

- Empty string -> clear parse error
- Just `;` -> error
- Just `()` with no semicolon -> error
- Unbalanced open: `(A,B` -> error with position
- Unbalanced close: `A,B)` -> error with position
- Missing semicolon: `(A,B)` -> error
- Double semicolon: `(A,B);;` -> parse first tree, ignore rest (or error -- decide)
- Binary data / NUL bytes -> error
- Deep nesting (10K levels) -> stack behavior documented (pest uses iterative parsing, should handle)

### Pathological -- performance

- Large tree: parse all `data/*/tree.nwk` dataset trees (up to 2000 leaves)
- Very large synthetic tree: 10K+ leaves (parsing time reasonable, no OOM)

### Golden master

- Parse every `data/*/tree.nwk` with both bio and pest parsers, compare resulting `Graph` structure: same topology, same names, same edge count
- Branch lengths: pest f64 values should be at least as precise as bio f32 values (compare against known f64 reference)

### Regression

- All existing tests calling `nwk_read`, `nwk_read_file`, `nwk_read_str` still pass
- All existing golden master tests still pass (branch length differences within f32->f64 improvement are expected and acceptable)

## Related issues

- Source: [kb/issues/M-io-bio-crate-newick-rejects-comments.md](../issues/M-io-bio-crate-newick-rejects-comments.md)
- Source: [kb/issues/N-nwk-branch-length-f32-precision-loss.md](../issues/N-nwk-branch-length-f32-precision-loss.md)
