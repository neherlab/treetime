# Normalize mathematical notation and symbols

Convert affected algorithm documentation and rustdoc to consistent rendered mathematics with complete symbol declarations.

## Required changes

- Replace ASCII/code-block equations with inline or display KaTeX.
- Add one `where` clause for every display equation’s nontrivial symbols.
- Define reroot distances, variances, indices, and summation domains.
- Update the canonical generator so generated reference output renders the same KaTeX equations and symbol declarations as its source documentation; remove the duplicated ASCII representation.

## Validation

- Run documentation generation and citation/math checks.
- Inspect rendered equations for valid KaTeX and unambiguous symbol scope.
- Regenerate checked-in CLI reference through its canonical generator.

## Related issues

- Source: [kb/issues/N-doc-mathematical-notation-inconsistent.md](../issues/N-doc-mathematical-notation-inconsistent.md)
