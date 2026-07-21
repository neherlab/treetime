# v0 Errata

Defects, oversights, and undocumented simplifications in v0 (Python) that v1 (Rust) does not reproduce. One file per erratum.

Distinct from:

- [_raw](../_raw/) - human-produced source material (specs, papers, notes), read-only for AI
- [algo](../algo/README.md) - algorithm documentation, scientific background, implementation status
- [decisions](../decisions/README.md) - deliberate v1 divergences from v0
- [features](../features/README.md) - feature parity checklist (done/partial/missing)
- [issues](../issues/README.md) - v1 defects, missing features, unintentional v0 divergences
- [proposals](../proposals/README.md) - undecided design documents with options and tradeoffs
- [reports](../reports/README.md) - research reports on algorithms, optimization methods, and implementation analysis
- [tickets](../tickets/README.md) - actionable implementation instructions derived from issues and proposals

## When to add an erratum

Create an entry when v1 diverges from v0 because v0 is wrong. Check this directory before matching v0 behavior.

Evidence standard (at least two):

- The function accepts a parameter for the correct behavior but the caller omits it
- Adjacent code (same module, same author) handles the case correctly
- The docstring or comments describe the correct behavior but code does not implement it
- The behavior contradicts the cited scientific reference
