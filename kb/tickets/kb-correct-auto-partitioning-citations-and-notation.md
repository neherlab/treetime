# Correct auto-partitioning citations and notation

Repair the citation graph, bibliographic metadata, and symbol declarations in [`kb/reports/auto-partitioning.md`](../reports/auto-partitioning.md) without changing its recommendations.

## Acceptance criteria

- Add citation anchors to the glossary citations for references 1, 3, 6, 8, and 9.
- Add first-use-linked glossary entries for CAT, AIC, and BIC; preserve the existing GTR and PMSF first-use links.
- Number repeated citation anchors monotonically and make every reference backlink resolve to its corresponding use.
- Use `et al.` for in-text citations with three or more authors.
- Format every reference title in sentence case.
- Correct the Liu article to Qin Liu, Michael A. Charleston, Shane A. Richards, and Barbara R. Holland; 2023; _Systematic Biology_ 72(1):92-105.
- Preserve the published diacritics in Hervé Philippe and Simon Tavaré.
- Declare $\pi$, $L$, and the Kronecker delta before use.
- Use distinct symbols for alphabet size and quantile-bin count and update every dependent equation consistently.
- Validate all internal links and DOI targets.

## Related issues

- Source: [kb/issues/N-kb-auto-partitioning-citation-metadata-and-notation.md](../issues/N-kb-auto-partitioning-citation-metadata-and-notation.md)
