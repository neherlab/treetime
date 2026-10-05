# Documentation references and source claims are not verifiable

Touched algorithm documents and proposals contain orphan references, missing anchors/backlinks, and external schema/implementation claims without authoritative sources.

## Remaining work

- Tree-format proposals and adapter rustdoc need pinned source links to the cloned Augur and UShER revisions, and to the PhyloXML specification, directly after each external behavior claim
- Every retained bibliography entry needs an inline citation; orphan entries without support are removed
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
