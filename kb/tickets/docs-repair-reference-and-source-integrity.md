# Repair reference and source integrity

Make every scientific reference and external behavior claim in the touched documentation verifiable.

## Required changes

- Replace dead/wrong Pearl, Cormen, and Brent links with verified authoritative records.
- Add inline citations for every retained bibliography entry; remove unsupported orphan entries.
- Add stable anchors/backlinks and renumber in first-citation order.
- Add pinned Augur, UShER, and PhyloXML source/spec links directly after external behavior claims.
- Replace truncated source references with full project paths or pinned permalinks.

## Validation

- Resolve every DOI/URL and every local code link.
- Check that each reference has at least one inline citation and each inline citation has a reference.
- Verify every pinned external path exists at the cited revision in `.repos/` using read-only Git commands.

## Related issues

- Source: [kb/issues/N-doc-reference-and-source-integrity.md](../issues/N-doc-reference-and-source-integrity.md)
