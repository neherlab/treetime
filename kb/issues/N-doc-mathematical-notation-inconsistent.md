# Mathematical notation is inconsistent and leaves symbols undeclared

Algorithm Markdown, proposals, generated CLI reference text, and reroot rustdoc use ASCII or code spans for mathematical expressions despite KaTeX support. Reroot sufficient-statistics formulas introduce symbols without defining their meaning or summation domain.

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

## Related issues

- [N-doc-reference-and-source-integrity.md](N-doc-reference-and-source-integrity.md)
