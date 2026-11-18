# Unimplemented algorithm inventory has incomplete citation structure

[`kb/algo/unimplemented.md`](../algo/unimplemented.md) has citation and mathematical-notation defects:

- The document has no glossary linking technical terms to first use
- The references ([kb/algo/unimplemented.md#L255-L259](../algo/unimplemented.md#L255-L259)) are bullets rather than an ordered list, their numbers skip `ref-2` (`ref-1`, `ref-3`, `ref-4`), and the article titles use title case instead of sentence case
- The GTR random-model step ([kb/algo/unimplemented.md#L169](../algo/unimplemented.md#L169)) opens the math span `$\tfrac12(W+W^\mathsf{T})` and closes it with a backtick, so the formula does not render

> [!IMPORTANT]
> **Investigation required.** Check whether the remaining math symbols ($N$ at L64, $f$, $g$, $\tau$, $t$ at L135, $\pi$, $W$, $\mu$ at L169 and L248) are declared at first use.

These defects make the cross-algorithm inventory harder to audit but do not imply that any unimplemented algorithm is ready for implementation.

## Fix

- Add a glossary before the references and link each technical-term entry to its first in-text use
- Convert the references to an ordered list in first-use order with consecutive numbers, with article titles in sentence case
- Close the math span at L169 and declare every symbol before first use
- Keep source-code links, feature status, and unresolved decision text unchanged

## Validation

- All internal anchors, repository links, and DOI targets resolve
