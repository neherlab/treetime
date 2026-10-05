# Mathematical notation is inconsistent and leaves symbols undeclared

Algorithm Markdown, proposals, generated CLI reference text, and reroot rustdoc use ASCII or code spans for mathematical expressions despite KaTeX support. Reroot sufficient-statistics formulas introduce symbols without defining their meaning or summation domain. For example, the skyline stiffness help in [docs/docs/reference.md#L173](../../docs/docs/reference.md#L173) writes `(stiffness/2) * Σ (ln(Tc_{i+1}/Tc_i))^2` in a code span.

## Scope

- `kb/algo/ancestral.md`
- `kb/algo/reroot.md`
- `kb/algo/timetree.md`
- `kb/proposals/reroot-generic-scoring-architecture.md`
- `docs/docs/reference.md`
- reroot and optimize rustdoc under `packages/treetime/src`

## Potential solutions

- O1. Convert authored Markdown/rustdoc directly and repair the reference generator for generated text.
- O2. Add a documentation transform after generation. This obscures the source notation and can diverge from authored documentation.

## Recommendation

Render Markdown and Rustdoc mathematics as KaTeX. Declare every nontrivial symbol once immediately before or after its first equation; keep Rust identifiers in code spans only when referring to code entities.

- Replace ASCII and code-block equations with inline or display KaTeX
- Give every display equation one `where` clause that declares its nontrivial symbols
- Define the reroot distances, variances, indices, and summation domains
- Update the canonical generator so the generated reference renders the same KaTeX equations and symbol declarations as its source documentation, without a duplicate ASCII form

## Validation

- Run documentation generation and the citation and math checks
- Inspect the rendered equations for valid KaTeX and unambiguous symbol scope
- Regenerate the checked-in CLI reference through its canonical generator

## Related issues

- [N-doc-reference-and-source-integrity.md](N-doc-reference-and-source-integrity.md)
