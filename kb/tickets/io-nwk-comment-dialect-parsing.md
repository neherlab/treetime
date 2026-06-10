# Add Newick comment dialect parsing

Extend the pest grammar from ticket #2 with BEAST, NHX, and plain comment productions. Inline per-comment dialect detection. Populate `BTreeMap<String, String>` on `NodeFromNwk::from_nwk`.

## Done (`util-newick` crate)

The `util-newick` crate (`packages/util-newick/`) implements the standalone dialect parser:

- Inline per-comment dialect detection: `&&NHX:` -> NHX, `&` -> BEAST, else -> plain (discarded or stored as raw comment)
- BEAST value types: scalars, quoted strings, `{...}` arrays, booleans (case-insensitive), bare-key flags
- Node vs branch annotation placement: BEAST2 canonical (before `:` = node, after `:` = branch), MrBayes (after length = branch)
- NHX colon-separated key=value pairs
- Typed `NewickValue` enum (`String`, `Number`, `Boolean`, `Array`) instead of flat `BTreeMap<String, String>`

## Remaining (integration into `treetime-io`)

- Wire `util-newick` parsed annotations into `treetime-io` node/edge `BTreeMap` population at `nwk.rs:73`
- Map `NewickValue` variants to `treetime-io`'s string-valued `BTreeMap<String, String>`
- Optional `--nwk-dialect` CLI override for malformed files

## Original scope

- Add comment productions to the pest grammar (~30 lines): Comment, CommentBody, BeastPairs, NhxPairs, Key, Value, Array, QuotedStr, Scalar, PlainText
- Inline dialect detection per comment: first chars after `[` determine mode (`&&NHX:` -> NHX, `&` -> BEAST, else -> discard)
- Populate `BTreeMap<String, String>` at `nwk.rs:73` (currently hardcoded empty)
- Handle value types: scalars, quoted strings, `{...}` arrays, `#RRGGBB` colors
- Handle node vs branch annotation placement (BEAST2: node meta before `:`, branch meta after `:`)
- Optional `--nwk-dialect` CLI override for malformed files

## Grammar sketch

```pest
Comment     = _{ "[" ~ CommentBody ~ "]" }
CommentBody = { "&&NHX:" ~ NhxPairs | "&" ~ BeastPairs | PlainText }
BeastPairs  = { BeastPair ~ ("," ~ BeastPair)* }
BeastPair   = { Key ~ "=" ~ Value }
NhxPairs    = { NhxPair ~ (":" ~ NhxPair)* }
NhxPair     = { Key ~ "=" ~ Value }
Value       = { Array | QuotedStr | Scalar }
Array       = { "{" ~ Scalar ~ ("," ~ Scalar)* ~ "}" }
PlainText   = { (!"[" ~ !"]" ~ ANY)* }
```

## Tests

### Happy paths -- BEAST dialect

- Single attribute: `[&mutations="A55G"]` -> `{"mutations": "A55G"}`
- Multiple attributes: `[&mutations="A55G",country="USA"]` -> both keys populated
- Numeric value: `[&rate=0.003]` -> `{"rate": "0.003"}`
- Boolean-like: `[&is_root=true]` -> `{"is_root": "true"}`
- Array value: `[&colors={#FF0000,#00FF00}]` -> `{"colors": "{#FF0000,#00FF00}"}` (stored as string)
- Quoted string with spaces: `[&label="New York"]` -> `{"label": "New York"}`
- Color value: `[&color=#FF0000]` -> `{"color": "#FF0000"}`
- Node meta before colon: `(A[&node_attr=1]:0.1,B:0.2);` -> annotation on node
- Branch meta after colon before length: `(A:0.1[&branch_attr=2],B:0.2);` -> annotation on branch/edge

### Happy paths -- NHX dialect

- Single attribute: `[&&NHX:S=human]` -> `{"S": "human"}`
- Multiple attributes: `[&&NHX:S=human:D=Y:B=100]` -> all three keys populated
- Standard tags: S (species), D (duplication), B (bootstrap), T (taxonomy ID), E (EC number)

### Happy paths -- plain comments

- Plain text: `[this is a comment]` -> discarded, empty BTreeMap
- Empty comment: `[]` -> discarded

### Happy paths -- dialect detection

- BEAST detection: `[&` prefix -> BEAST parsing
- NHX detection: `[&&NHX:` prefix -> NHX parsing
- Plain detection: `[` without `&` -> discard
- Mixed dialects in one file: `(A[&k=v]:0.1,B[&&NHX:S=human]:0.2);` -> each comment parsed independently

### Edge cases -- values

- BEAST value with comma inside quotes: `[&label="A,B"]` -> `{"label": "A,B"}`
- BEAST value with equals inside quotes: `[&formula="x=y"]` -> `{"formula": "x=y"}`
- NHX value with embedded colon: ambiguous -- document behavior (truncate at colon or include?)
- Unquoted BEAST value with path chars: `[&path=/usr/local]` -> `{"path": "/usr/local"}`
- Empty value: `[&key=]` -> `{"key": ""}`
- Value with leading/trailing whitespace: `[&key= value ]` -> trimmed or preserved? Document
- Key with dots: `[&posterior.prob=0.95]` -> valid key
- Key with hyphens: `[&node-type=internal]` -> valid key
- Key with underscores: `[&num_date=2020.5]` -> valid key

### Edge cases -- structure

- Multiple comments on one node: `(A[&k1=v1][&k2=v2]:0.1)` -> merged BTreeMap
- Comment on root: `(A,B)[&root=true];` -> root node annotated
- Comment on leaf with no length: `(A[&k=v],B);` -> annotation without branch length
- 50+ attributes in one comment
- Annotation key collision between NodeToNwk::nwk_comments() and CommentProviders: later source wins (verify merge order at nwk.rs:239)

### Edge cases -- malformed

- `[&` with no closing `]` -> parse error with position
- `[&&NHX` with no colon -> parse error or treat as plain comment?
- `[&key]` with no `=value` -> parse error or ignore?
- `[&=value]` with no key -> parse error or ignore?
- `[&key=value` missing closing bracket -> parse error with position
- Nested brackets: `[&label="a[b]c"]` -> behavior documented

### Pathological

- Multi-megabyte comment content -> parse without OOM
- Binary data inside comment brackets -> error, not crash
- 1000+ annotations on a single node -> performance acceptable

### Round-trip

- Parse BEAST-annotated tree, write back with `NwkStyle::Annotated`, compare: parsed keys/values survive round-trip
- Parse NHX-annotated tree, write back with `NwkStyle::Nhx`, compare
- Parse plain-comment tree, write back: comments absent from output (expected)

### Golden master

- Parse annotated trees from external tools (FigTree, BEAST output, IQ-TREE output) if sample files available in test fixtures
- Parse v1's own annotated output (from pre-ticket-#1 behavior) -> verify annotations recovered

### Regression

- Trees without comments: identical parse result to ticket #2 (no regression from adding comment grammar)

## Related issues

- Source: [kb/issues/M-io-bio-crate-newick-rejects-comments.md](../issues/M-io-bio-crate-newick-rejects-comments.md)
- Dialect grammars: [kb/reports/newick-annotation-dialects.md](../reports/newick-annotation-dialects.md)
