# Feature Inventory - TreeTime v0/v1

> **Status**: reviewed against source code
>
> **Valid for commit**: `404d8b77` on 2026-03-10
>
> **Scope**: v0/v1 feature parity tracking. Combines domain taxonomy with CLI wiring status.

## Legend

Uses [Obsidian checkbox statuses](https://publish.obsidian.md/tasks/Getting+Started/Statuses):

- `[x]` implemented in v1 with v0 parity (or v1-only feature)
- `[/]` partial implementation, stubbed, or only partly wired
- `[ ]` missing, parsed but not wired, or documented but unimplemented

## Command Map

| Command     | Status | Notes                                      |
| ----------- | ------ | ------------------------------------------ |
| `ancestral` | [x]    | Parsimony and marginal reconstruction      |
| `clock`     | [x]    | Regression and rerooting                   |
| `timetree`  | [x]    | Full inference pipeline                    |
| `optimize`  | [x]    | v1-only branch-length optimization         |
| `prune`     | [x]    | v1-only tree pruning                       |
| `mugration` | [x]    | Marginal reconstruction with iterative GTR |
| `homoplasy` | [ ]    | Unimplemented                              |

## Domain Pages

| Domain                        | Page                                   | v0 Features | v1 Done | v1 Missing | v1 Only |
| ----------------------------- | -------------------------------------- | ----------- | ------- | ---------- | ------- |
| Ancestral Reconstruction      | [ancestral.md](ancestral.md)           | 16          | 10      | 6          | 0       |
| Clock Inference               | [clock.md](clock.md)                   | 17          | 14      | 3          | 0       |
| Timetree Inference            | [timetree.md](timetree.md)             | 30          | 21      | 9          | 0       |
| Homoplasy Analysis            | [homoplasy.md](homoplasy.md)           | 9           | 0       | 9          | 0       |
| Mugration                     | [mugration.md](mugration.md)           | 10          | 10      | 0          | 0       |
| ARG Inference                 | [arg.md](arg.md)                       | 4           | 0       | 4          | 0       |
| Branch Length Optimization    | [optimize.md](optimize.md)             | 5           | 3       | 5          | 2       |
| Pruning                       | [prune.md](prune.md)                   | 0           | 0       | 0          | 7       |
| GTR Substitution Models       | [gtr.md](gtr.md)                       | 12          | 7       | 5          | 0       |
| Alphabet System               | [alphabet.md](alphabet.md)             | 6           | 5       | 1          | 0       |
| Probability Distributions     | [distributions.md](distributions.md)   | 15          | 5       | 10         | 0       |
| Representation and Partitions | [representation.md](representation.md) | 12          | 12      | 0          | 0       |
| Sequence Primitives           | [primitives.md](primitives.md)         | 7           | 7       | 0          | 0       |
| I/O Formats                   | [io.md](io.md)                         | 11+13       | 9+6     | 2+7        | 0+4     |
| Date Parsing                  | [dates.md](dates.md)                   | 5           | 4       | 1          | 0       |

## Cross-Command Notes

- `ancestral`, `timetree`, `homoplasy` expose `MethodAncestral`
- `clock` and `timetree` share the same clock-regression and reroot machinery
- `timetree` reuses `ancestral::marginal` and `clock::*` internals rather than maintaining a separate reconstruction stack
- `optimize` and `prune` are v1-only commands
- Several help strings describe planned behavior not wired in current command code

## Cross-References

- [Algorithm Inventory](../algo/README.md) - implementation details
- [Decisions](../decisions/README.md) - deliberate design choices
- [Issues](../issues/README.md) - bugs and missing features
- [Test Inventory](../tests/README.md) - test coverage
