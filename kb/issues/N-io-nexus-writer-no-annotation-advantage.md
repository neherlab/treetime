# Nexus writer provides no annotation advantage over Newick

The Nexus writer `nex_write_with` at [packages/treetime-io/src/nex.rs#L91-L137](../../packages/treetime-io/src/nex.rs#L91-L137) wraps the same Newick string (produced by `nwk_write_with`) in TAXA/TREES blocks. It does not use any NEXUS-specific mechanism for annotations.

The NEXUS specification defines no per-node or per-branch annotation mechanism. Comments `[...]` are formally discardable per the spec. The `[&key=value]` annotations in our `.nexus` output are BEAST-convention comments inside the embedded Newick string, not NEXUS features.

The `.nexus` and `.nwk` outputs carry identical annotation content. Prior design discussions incorrectly assumed `.nexus` "carries annotations" as a format advantage.

## Current Nexus output structure

```nexus
#NEXUS
Begin Taxa;
  Dimensions NTax=N;
  TaxLabels name1 name2 ...;
End;
Begin Trees;
  Tree tree1=(A:0.1[&mutations="G42T"],B:0.2)root;
End;
```

Missing NEXUS features that other tools use:

- No `TRANSLATE` table (taxon number-to-name mapping, reduces file size for large trees)
- No `PROPERTIES` command (e.g. `rooted=yes`)
- No FigTree-style custom blocks for display settings

## Related

- [M-io-nwk-writer-annotation-format-nonstandard.md](M-io-nwk-writer-annotation-format-nonstandard.md)
- [kb/reports/newick-annotation-dialects.md](../reports/newick-annotation-dialects.md)
