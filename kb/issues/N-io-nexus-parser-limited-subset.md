# Nexus parser handles a limited subset

The `util-newick` crate Nexus parser covers the TREES block subset needed for phylogenetic tree I/O. It does not implement a full Nexus grammar parser.

## What is supported

- `Begin Trees;` / `End;` block extraction (case-insensitive)
- `Translate` integer-to-name tables with quote-aware comma splitting
- `Tree name = [&R|U]? newick;` entries with quoted name handling
- Multiple trees per TREES block
- Unknown blocks silently skipped (FigTree, Data, Characters, etc.)

## What is not supported

- `Begin Characters;`, `Begin Data;`, `Begin Sets;`, `Begin Assumptions;` -- content ignored
- Nexus comments `[...]` at the block level (only Newick-embedded comments are handled)
- `UTREE`, `TREE * name =` syntax variants
- Semicolons inside double-quoted labels in tree commands
- Malformed blocks produce empty results rather than parse errors

## Why this is acceptable

Nexus is a 1997 format with no formal revision. The only tree-relevant content is inside TREES blocks. A full Nexus grammar covering CHARACTERS, ASSUMPTIONS, SETS, etc. would be a separate crate with no benefit for tree I/O.

## What would change this

If the project needs to read Nexus data matrices or other non-tree blocks, a formal Nexus parser (pest grammar or similar) should replace the string scanner in `nexus.rs`.
