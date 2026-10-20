# Documentation references and source claims are not verifiable

Algorithm documents and proposals in `kb/` contain orphan references, references without backlinks, and claims about external formats without authoritative sources.

## Remaining work

- [kb/proposals/tree-format-model-embedding.md#L11-L12](../proposals/tree-format-model-embedding.md#L11-L12) states that the Auspice v2 schema and the UShER MAT protobuf have no model field and that augur keeps the model out of the Auspice JSON, without a source link. Each claim needs a pinned link to the cloned Augur, Auspice, or UShER revision directly after it
- Reference entries without an inline citation (orphans): `kb/algo/ancestral.md` (`ref-2`, `ref-3`, `ref-6`), `kb/algo/mugration.md` (`ref-0`, `ref-5`, `ref-6`), `kb/algo/timetree.md` (`ref-8`, `ref-13`), `kb/proposals/config-file-multi-partition.md` (`ref-1` to `ref-7`). Each needs an inline citation where the document relies on it, or removal
- Reference entries without a backlink to their citation: the orphans above, and `ref-1` to `ref-7` of `kb/proposals/graph-pure-transform-architecture.md`
- References sections without anchors for inline citations: `kb/proposals/mugration-full-reconstruction-per-iteration.md`, `kb/proposals/optimize-convergence-and-robustness.md`, `kb/proposals/optimize-indel-rate-per-iteration-reestimation.md`
- References need stable anchors and backlinks, numbered in first-citation order
- Truncated source references need full project paths or pinned permalinks

The Pearl, Cormen, and Brent records in `kb/algo/` already point to verified sources: the Pearl DOI `10.1016/C2009-0-27609-4` (`kb/algo/timetree.md`, `kb/algo/ancestral.md`), the MIT Press fourth-edition page for Cormen (`kb/algo/graph.md`), and the author-maintained bibliography for Brent (`kb/algo/clock.md`, `kb/algo/optimization.md`).

## Potential solutions

- O1. Resolve authoritative metadata and normalize citations/source links in one coordinated edit.
- O2. Repair only broken links. This leaves orphan references and unsupported external claims unverifiable.

## Recommendation

Repair metadata against authoritative sources, then normalize inline citations, reference anchors/backlinks, citation order, and full-path code references in one pass. Pinned external code links must include the exact inspected commit.

## Validation

- Resolve every DOI, URL, and local code link
- Check that each reference has at least one inline citation and each inline citation has a reference
- Verify that every pinned external path exists at the cited revision in `.repos/`, using read-only Git commands

## Related issues

- [N-doc-mathematical-notation-inconsistent.md](N-doc-mathematical-notation-inconsistent.md)
