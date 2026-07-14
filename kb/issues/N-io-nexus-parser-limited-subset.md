# Nexus parser handles a limited subset

> [!IMPORTANT]
> **Scope decision required.** The current tree-only subset is implemented deliberately, but its acceptance as the complete Nexus contract has not been approved. Keep this as an issue until the supported syntax and malformed-input behavior are explicitly decided.

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

## Current scope rationale

TreeTime currently needs tree content from `TREES` blocks. Supporting unrelated Nexus data matrices would expand the I/O boundary substantially, while silent acceptance of malformed tree blocks and unsupported tree-command variants remains observable within the current scope.

## Decision axes

- **Documented subset vs broader grammar:** retain a tree-only parser or adopt a general Nexus grammar.
- **Unsupported syntax:** reject unsupported tree syntax explicitly or continue skipping it.
- **Malformed input:** return actionable parse errors or preserve empty-result behavior.
- **Non-tree blocks:** ignore them as out of scope or represent them in a broader data model.

No implementation ticket is ready until these axes are decided.
