# Repair missing knowledge-base link targets

Correct every missing internal Markdown target listed in the source issue.

## Acceptance criteria

- Trace moved source files with Git history and update links to their current tracked locations and relevant line anchors.
- Replace links to deleted issues, tickets, and proposals with their live successor when one exists; otherwise remove the stale statement or rewrite it to stand alone.
- Replace links to untracked build artifacts with stable reproduction instructions or tracked evidence.
- Preserve the scientific and design meaning of each surrounding paragraph.
- Run a repository-wide Markdown target validator and confirm that every tracked relative target exists.
- Do not modify `kb/_raw/`.

## Related issues

- Source: [kb/issues/N-kb-internal-links-have-missing-targets.md](../issues/N-kb-internal-links-have-missing-targets.md)
