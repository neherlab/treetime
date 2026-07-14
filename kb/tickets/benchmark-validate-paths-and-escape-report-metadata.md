# Validate benchmark paths and escape report metadata

Separate benchmark command, filesystem, and Markdown data boundaries.

## Required changes

- Remove output-path interpolation from shell-valued Hyperfine arguments.
- Parse dataset identifiers without absolute paths, parent traversal, or separators.
- Normalize output paths and prove every generated path remains under the configured root.
- Escape metadata used in Markdown text, tables, and links.

## Validation

- Shell metacharacter, whitespace, absolute path, traversal, separator, and Markdown injection cases.
- Whole output-tree comparison proving no file escapes the root.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/H-benchmark-paths-cross-shell-filesystem-and-markdown-boundaries.md](../issues/H-benchmark-paths-cross-shell-filesystem-and-markdown-boundaries.md)
