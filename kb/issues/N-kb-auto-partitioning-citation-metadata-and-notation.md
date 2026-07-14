# Auto-partitioning report has citation metadata and notation defects

[`kb/reports/auto-partitioning.md`](../reports/auto-partitioning.md) contains independently checkable documentation defects that impede source verification:

## Citation graph

- Glossary citations for references 1, 3, 6, 8, and 9 lack their own in-text citation anchors.
- CAT, AIC, and BIC lack glossary entries linked to first use; the existing GTR and PMSF entries already have first-use links.
- Repeated citation anchors for reference 3 are ordered `cite-3`, `cite-3c`, `cite-3b`; repeated anchors for reference 15 start unsuffixed and then use `b`.
- References 1, 3, 6, 8, and 9 lack glossary backlinks.
- Three-or-more-author prose citations do not consistently use `et al.`.
- Reference titles do not consistently use sentence case.

## Metadata

- The Liu paper is a 2023 _Systematic Biology_ article, volume 72 issue 1, pages 92-105, by Qin Liu, Michael A. Charleston, Shane A. Richards, and Barbara R. Holland.
- Hervé Philippe's and Simon Tavaré's names require their published diacritics.

## Mathematical notation

- $\pi$ is used before its equilibrium-frequency meaning is declared.
- $L$ and the Kronecker delta are not declared.
- One symbol is reused for alphabet size and quantile-bin count.

These are documentation defects; they do not change a scientific or implementation decision.

## Related tickets

- [kb/tickets/kb-correct-auto-partitioning-citations-and-notation.md](../tickets/kb-correct-auto-partitioning-citations-and-notation.md)
