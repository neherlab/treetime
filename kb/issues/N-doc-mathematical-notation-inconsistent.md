# Mathematical notation is inconsistent and leaves symbols undeclared

Algorithm Markdown and the generated CLI reference write mathematical expressions as ASCII in code spans and code blocks, although the Markdown documents render KaTeX. Symbols in these expressions are often not declared. Examples:

- [kb/algo/ancestral.md#L48-L52](../algo/ancestral.md#L48-L52): the partial likelihood `w_k(X) = [sum_Y P(X->Y|t_i) * w_i(Y)] * ...` in a code block
- [kb/algo/timetree.md#L90-L92](../algo/timetree.md#L90-L92): the skyline cost in a code block; the confidence interval formulas at [kb/algo/timetree.md#L252-L257](../algo/timetree.md#L252-L257) and [kb/algo/timetree.md#L267-L270](../algo/timetree.md#L267-L270) as well
- [docs/docs/reference.md#L187](../../docs/docs/reference.md#L187): the skyline stiffness help writes `(stiffness/2) * Σ (ln(Tc_{i+1}/Tc_i))^2` in a code span. The text comes from the doc comment of the `--skyline-stiffness` flag [packages/app-commands/src/commands/timetree/args.rs#L276](../../packages/app-commands/src/commands/timetree/args.rs#L276), which `--help` also prints in a terminal

## Scope

- `kb/algo/ancestral.md`
- `kb/algo/timetree.md`
- The CLI help of `--skyline-stiffness` and the generated `docs/docs/reference.md`

> [!IMPORTANT]
> **Decision required.** The CLI reference is generated from the clap doc comments, which `--help` prints as plain text in a terminal. KaTeX in these comments renders in `docs/docs/reference.md` but shows as raw TeX (`$\frac{s}{2} \sum_i ...$`) in `--help`. Options:
>
> - **KaTeX in the doc comments**: the reference renders the formula; `--help` shows raw TeX
> - **Plain text in the doc comments**: `--help` stays readable; the reference keeps a plain-text formula
> - **Plain text in the doc comments and KaTeX in the reference**: the reference generator replaces the formula, so the help text has two sources

## Fix

For the KaTeX documents in `kb/algo/`:

- Replace ASCII and code-block equations with inline or display KaTeX
- Give every display equation one `where` clause that declares its nontrivial symbols
- Keep Rust identifiers in code spans only when the text refers to code entities

For the CLI help, apply the option that the decision above selects.

## Validation

- Inspect the rendered equations for valid KaTeX and unambiguous symbol scope
- Regenerate the checked-in CLI reference with `just gen`, and check `treetime timetree --help` in a terminal

## Related issues

- [N-doc-reference-and-source-integrity.md](N-doc-reference-and-source-integrity.md)
