# Nexus output quotes taxon labels by Newick rules

`fn nex_write()` in [packages/treetime-io/src/nex.rs](../../packages/treetime-io/src/nex.rs) writes the `TaxLabels` of the `Taxa` block with `fn write_label()` of `util-newick`, which quotes only the Newick punctuation. The Nexus standard (Maddison et al. 1997) also treats `{ } / \ = * " ` + - < >` as punctuation, so a label such as `hCoV-19/USA/CA-1/2020` is written unquoted. `NTax` counts every named leaf, so a tree with repeated leaf names lists the same label more than once.

> [!IMPORTANT]
> **Investigation required.** Biopython 1.88 reads both forms. No reader that rejects them has been found yet; check FigTree, PAUP* and DendroPy before changing the output, because the change alters every `.nexus` output of datasets with such names.

The Nexus writer of `util-newick` (`fn write_nexus_word()` in [packages/util-newick/src/nexus.rs](../../packages/util-newick/src/nexus.rs)) already quotes by the Nexus rules.
