# Unimplemented algorithm inventory has incomplete citation structure

[`kb/algo/unimplemented.md`](../algo/unimplemented.md) has citation and mathematical-notation defects across multiple algorithm entries:

- The document has no glossary linking technical terms to first use.
- Reference entries are formatted as bullets rather than an ordered reference list, and the article titles use title case instead of sentence case ([kb/algo/unimplemented.md#L321-L328](../algo/unimplemented.md#L321-L328)).
- The per-site rate variation entry uses $L_a$, $w_k$, $r_k$, $k$, $Q$, $V$, $\lambda$, and $t$ without declarations at first use ([kb/algo/unimplemented.md#L93](../algo/unimplemented.md#L93), [#L104](../algo/unimplemented.md#L104)). The same entry writes the per-site rate as both $\mu^a$ and $\mu_a$.

These defects make the cross-algorithm inventory harder to audit but do not imply that any unimplemented algorithm is ready for implementation.

## Fix

- Add a glossary before the references and link each technical-term entry to its first in-text use
- Convert the references to an ordered list in first-use order, with article titles in sentence case
- Declare every symbol before first use, with distinct symbols for distinct concepts and one notation for the per-site rate
- Keep source-code links, feature status, and unresolved decision text unchanged

## Validation

- All internal anchors, repository links, and DOI targets resolve
