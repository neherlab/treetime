# Knowledge Base

Shared knowledge base (KB). AI agents and humans collaborate here: documenting progress, tracking issues, recording design decisions, and maintaining project knowledge.

> 💡 ## LLM wiki pattern
>
> This knowledge base follows the [LLM wiki](https://gist.github.com/karpathy/442a6bf555914893e9891c11519de94f) pattern: raw source material is compiled by AI agents into structured, cross-linked knowledge articles. The organizing structure is human-defined. AI maintains content within that structure.

## Motivation

- AI continuity. Sessions start with zero memory. KB persists decisions, defects, and findings so knowledge compounds instead of being rediscovered.
- Human onboarding. Synthesized view of algorithms, design rationale, open problems, and research context.
- Multi-agent coordination. Parallel agents working on different tickets share state through issues, decisions, and proposals.
- Scientific rigor. Math-heavy phylogenetics code requires tracking what is correct, what diverges from theory, and what diverges from the reference implementation.
- Preventing rework. Known issues and errata prevent re-investigating solved problems or re-introducing fixed bugs.
- Decision traceability. Design choices are documented with rationale, not implicit in code.

## Directories

| Directory                  | Description                                                                                                                                                                                                                                               |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`_raw/`](_raw/)           | Human-produced source material (specifications, papers, notes). Read-only for AI.                                                                                                                                                                         |
| [`algo/`](algo/)           | Algorithm documentation: scientific background, implementation status, v0/v1 locations                                                                                                                                                                    |
| [`decisions/`](decisions/) | Deliberate v1 design choices with rationale (one file per decision)                                                                                                                                                                                       |
| [`features/`](features/)   | Feature parity checklist: `[x]` done, `[/]` partial, `[ ]` not done                                                                                                                                                                                       |
| [`issues/`](issues/)       | Concrete problems. Severity-prefixed (H/M/N). The working list agents consult before domain work. PREFER independent issues, but entangled problems may share a file when splitting would lose clarity                                                    |
| [`proposals/`](proposals/) | Undecided design documents analyzing a problem space with options and tradeoffs. Source material for issues -- every actionable item in a proposal must be extracted into a separate issue so it is not lost when the proposal is no longer actively read |
| [`reports/`](reports/)     | Research reports on algorithms, optimization methods, and implementation analysis                                                                                                                                                                         |
| [`tickets/`](tickets/)     | Implementation instructions for a coding agent. One task per file, executable in one session without further research or decisions. Derived from issues and proposals when the implementation path is fully decided                                       |
| [`v0-errata/`](v0-errata/) | Defects in v0 that v1 correctly avoids (2+ evidence sources required)                                                                                                                                                                                     |

## Structure

[`_raw/`](_raw/) contains source material. All other directories contain AI-maintained derived knowledge. Source code is ground truth. KB entries are guides, not substitutes for code verification.

When new material is added to [`_raw/`](_raw/), dependent articles in other directories should be reviewed and updated.

## Taxonomy

Every work item falls into exactly one category:

| Category                        | Directory                                                                     | Scope                                                                |
| ------------------------------- | ----------------------------------------------------------------------------- | -------------------------------------------------------------------- |
| Implemented, same behavior      | [`algo/`](algo/), [`features/`](features/)                                    | Feature with equivalent behavior to v0 or design spec                |
| Implemented, different behavior | [`decisions/`](decisions/)                                                    | Feature with deliberate divergence from v0, or intentionally removed |
| Not yet done                    | [`issues/`](issues/), [`algo/unimplemented.md`](algo/unimplemented.md)        | Bugs, missing features, stubs, dead flags, behavioral differences    |
| New in v1                       | [`proposals/`](proposals/) (pre-impl), [`decisions/`](decisions/) (post-impl) | Feature not in v0 or design specs                                    |
| v0 defective, v1 correct        | [`v0-errata/`](v0-errata/)                                                    | v0 defect that v1 does not reproduce                                 |

### Decision rules

- Feature working with same results: [`features/`](features/) `[x]`, algorithm in domain file
- Feature working with different results: [`decisions/`](decisions/) (one file with rationale)
- Feature intentionally removed: [`decisions/`](decisions/) (one file with rationale)
- Feature missing, stubbed, or broken: [`issues/`](issues/) (severity-prefixed file). If algorithmic, also in [`algo/unimplemented.md`](algo/unimplemented.md)
- Design-spec feature not yet implemented: [`issues/`](issues/) (severity-prefixed file)
- New v1 feature: [`proposals/`](proposals/) pre-implementation, then [`decisions/`](decisions/) post-implementation
- v0 behavior wrong, v1 correct: [`v0-errata/`](v0-errata/) (one file with evidence)

### Severity (issues only)

| Prefix | Severity   | Criteria                                                                        |
| ------ | ---------- | ------------------------------------------------------------------------------- |
| `H-`   | High       | Crashes, panics, blocks correct results, or missing standard expected feature   |
| `M-`   | Medium     | Wrong results under specific conditions, or missing feature affecting workflows |
| `N-`   | Negligible | Edge cases, niche missing features, weak assertions, cosmetic                   |

### Proposal lifecycle

- Proposal created during research session with ecosystem survey, design axes, options, tradeoffs
- User decides per axis. Decided items become tickets (if immediately implementable) or stay as issues (if further research needed)

### Ticket lifecycle

- Issue exists in [`issues/`](issues/). Ticket created in [`tickets/`](tickets/) with `## Related issues` linking back
- Ticket readiness: all design decisions made, implementation path clear, no open questions requiring user input. An issue with undecided design axes is not ready for a ticket
- Ticket executed. Both deleted if fully resolved
- Partial resolution: update both to reflect remaining work
