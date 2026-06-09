# Implement format adapter traits on command payload types

The format adapter traits (`PhyloxmlFromGraph`, `AuspiceWrite`, `UsherWrite`) are defined in treetime-io with zero implementations on analysis command payload types. Implement them so the output selection system (ticket #4) can dispatch to these formats.

## Scope

- `PhyloxmlFromGraph` on `NodeAncestral`/`EdgeAncestral`, `NodeTimetree`/`EdgeTimetree`, `NodeClock`/`EdgeClock`
- `AuspiceWrite` on `NodeAncestral`/`EdgeAncestral` (mutations), mugration payloads (discrete traits as colorings). `NodeTimetree` already has `TimetreeAuspiceWriter` at `timetree/output/auspice.rs`
- `UsherWrite` on `NodeAncestral`/`EdgeAncestral` and `NodeTimetree`/`EdgeTimetree` (require per-branch mutation data). Not meaningful for clock, prune, mugration

## Per-command availability after implementation

| Format          | ancestral | optimize | timetree | mugration | clock | prune |
| --------------- | --------- | -------- | -------- | --------- | ----- | ----- |
| auspice         | new       | no       | existing | new       | no    | no    |
| phyloxml        | new       | new      | new      | new       | new   | new   |
| mat-pb/mat-json | new       | new      | new      | no        | no    | no    |

## Locations

- PhyloXML traits: `packages/treetime-io/src/phyloxml.rs:195,257`
- Auspice traits: `packages/treetime-io/src/auspice.rs`
- UShER traits: `packages/treetime-io/src/usher_mat.rs`
- Existing auspice impl: `packages/treetime/src/commands/timetree/output/auspice.rs:46`

## Tests

### PhyloXML -- happy paths

- Write tree with names, branch lengths -> valid PhyloXML (parseable by quick-xml)
- Leaf names map to `<name>` element in `<clade>`
- Branch lengths map to `<branch_length>` element
- Mutations map to `<property>` elements (ancestral, timetree, optimize)
- Dates map to `<date>` element (timetree)
- Discrete traits map to `<property>` elements (mugration)
- Round-trip: write PhyloXML, read back via `phyloxml_to_graph`, compare topology and branch lengths

### PhyloXML -- edge cases

- Names with XML-special chars (`<`, `>`, `&`, `"`, `'`) -> properly escaped
- Unicode names (umlaut, CJK)
- Empty tree (root only, no children)
- Node with no name (unnamed internal)
- Node with no branch length

### Auspice (non-timetree) -- happy paths

- Ancestral: mutations in `branch_attrs.mutations`, divergence in `node_attrs.div`
- Mugration: trait values as `node_attrs.{attribute}` with `value` field, colorings populated in meta
- Output validates against Auspice v2 JSON schema
- Round-trip: write Auspice JSON, read back via `auspice_to_graph`, compare topology

### Auspice (non-timetree) -- edge cases

- Node with no mutations -> empty `branch_attrs.mutations`
- Mugration with unknown/missing trait value (`?`) -> handled in node_attrs
- Trait values with special JSON chars (quotes, backslashes) -> properly escaped
- 100+ mutations per branch -> all serialized

### UShER MAT -- happy paths

- Ancestral: per-branch mutations encoded in protobuf `mut` messages
- Mutation encoding: position, ref_nuc, par_nuc, mut_nuc (0=A, 1=C, 2=G, 3=T)
- Condensed nodes: identical sequences grouped (if applicable)
- Round-trip: write MAT protobuf, read back via `usher_to_graph`, compare mutations

### UShER MAT -- edge cases

- Branch with no mutations -> empty mutation list
- Ambiguous nucleotides in mutations (N, deletion `-`) -> encoding documented
- Indels -> not representable in MAT format, skipped or error? Document
- Tree with 10K+ nodes -> serialization performance acceptable

### Pathological (all formats)

- Empty graph (no nodes) -> error or valid empty document? Document per format
- Single-node graph (root only) -> valid minimal document
- Graph with 100K+ nodes -> no OOM, reasonable time

### Regression

- Existing `TimetreeAuspiceWriter` output unchanged
- Existing tests for timetree auspice JSON still pass

## Related issues

- Source: [kb/issues/N-io-write-graph-files-missing-formats.md](../issues/N-io-write-graph-files-missing-formats.md)
