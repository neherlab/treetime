# Known Issues

Tracks unintentional bugs, missing features, todos, behavioral differences from v0.

Distinct from:

- [_raw](../_raw/) - human-produced source material (specs, papers, notes), read-only for AI
- [algo](../algo/README.md) - algorithm documentation, scientific background, implementation status
- [decisions](../decisions/README.md) - deliberate v1 divergences from v0
- [features](../features/README.md) - feature parity checklist (done/partial/missing)
- [proposals](../proposals/README.md) - undecided design documents with options and tradeoffs
- [reports](../reports/README.md) - research reports on algorithms, optimization methods, and implementation analysis
- [tickets](../tickets/README.md) - actionable implementation instructions derived from issues and proposals
- [v0-errata](../v0-errata/README.md) - v0 defects that v1 correctly avoids

## Filename convention

- Files are prefixed with a severity letter so that letters H < M < N sort alphabetically in severity order:

| Prefix | Severity   | Criteria                                                                        |
| ------ | ---------- | ------------------------------------------------------------------------------- |
| `H-`   | High       | Crashes, panics, blocks correct results, or missing standard expected feature   |
| `M-`   | Medium     | Wrong results under specific conditions, or missing feature affecting workflows |
| `N-`   | Negligible | Edge cases, niche missing features, weak assertions, cosmetic                   |

- Domain prefix (`ancestral-`, `clock-`, `timetree-`, etc.) follows the severity letter, grouping related issues within each tier.
- Severity for missing features: a missing feature that most users expect (e.g., a standard phylogenetic capability) is High. A missing feature that affects specific workflows is Medium. A niche or rarely used missing feature is Negligible.

Issue name in the summary table must match the H1 heading in the linked file exactly.
