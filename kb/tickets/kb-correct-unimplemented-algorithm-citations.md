# Correct unimplemented algorithm inventory citations

Repair the citation graph and symbol declarations in [`kb/algo/unimplemented.md`](../algo/unimplemented.md) without changing feature status or deciding any open implementation contract.

## Acceptance criteria

- Add a glossary before the references and connect each technical-term entry to its first in-text use.
- Convert references to an ordered list in first-use order and use sentence case for article titles.
- Add in-text uses and backlinks for the Dempster and Felsenstein references where they support the implemented iterative-GTR entry, or remove references that have no supported use.
- Declare $L_a$, $w_k$, $r_k$, $k$, $Q$, $V$, $\lambda$, and $t$ before first use, with distinct symbols for distinct concepts.
- Preserve source-code links, feature status, and unresolved decision language.
- Validate all internal anchors, repository links, and DOI targets.

## Related issues

- Source: [kb/issues/N-kb-unimplemented-algorithm-citation-structure.md](../issues/N-kb-unimplemented-algorithm-citation-structure.md)
