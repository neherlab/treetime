# Documentation references and source claims are not verifiable

Touched algorithm documents and proposals contain orphan references, missing anchors/backlinks, dead or incorrect DOI links, and external schema/implementation claims without authoritative sources.

## Verified metadata failures

- Pearl’s book DOI should be `10.1016/C2009-0-27609-4`; the current DOI does not resolve.
- The cited Cormen DOI does not resolve; use the official MIT Press fourth-edition page.
- The Brent DOI resolves to an unrelated book; use the author-maintained bibliography.
- Tree-format proposals and adapter rustdoc need pinned source links to the cloned Augur and UShER revisions.

## Potential solutions

- O1. Resolve authoritative metadata and normalize citations/source links in one coordinated edit.
- O2. Repair only dead links. This leaves orphan references and unsupported external claims unverifiable.

## Recommendation

Repair metadata against authoritative sources, then normalize inline citations, reference anchors/backlinks, citation order, and full-path code references in one pass. Pinned external code links must include the exact inspected commit.

## Related issues

- [N-code-quality-conventions.md](N-code-quality-conventions.md)
