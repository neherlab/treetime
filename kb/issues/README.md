# Known Issues

Tracks unintentional bugs, missing features, todos, behavioral differences from v0.

Distinct from:

- [_raw](../_raw/) - human-produced source material (specs, papers, notes), read-only for AI
- [algo](../algo/README.md) - algorithm documentation, scientific background, implementation status
- [decisions](../decisions/README.md) - deliberate v1 divergences from v0
- [features](../features/README.md) - feature parity checklist (done/partial/missing)
- [proposals](../proposals/README.md) - undecided design documents with options and tradeoffs
- [reports](../reports/README.md) - research reports on algorithms, optimization methods, and implementation analysis
- [tickets](../tickets/README.md) - actionable implementation instructions derived from decided, implementation-ready issues
- [v0-errata](../v0-errata/README.md) - v0 defects that v1 correctly avoids

## Filename convention

- Files are prefixed with a severity letter so that letters H < M < N sort alphabetically in severity order:

| Prefix | Severity   | Criteria                                                                        |
| ------ | ---------- | ------------------------------------------------------------------------------- |
| `H-`   | High       | Crashes, data loss, incorrect scientific results, or blocked required behavior  |
| `M-`   | Medium     | Incorrect behavior under bounded conditions or a specified capability gap       |
| `N-`   | Negligible | Documentation, test, maintainability, or presentation defect with no demonstrated runtime effect |

- Domain prefix (`ancestral-`, `clock-`, `timetree-`, etc.) follows the severity letter, grouping related issues within each tier.
- Derive severity from specification language, user requirements, external evidence, or demonstrated impact. If the evidence does not distinguish a severity, do not infer one from assumed usage frequency.
