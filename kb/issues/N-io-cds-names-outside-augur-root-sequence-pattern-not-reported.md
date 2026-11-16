# CDS names outside the augur root-sequence pattern are not reported

The augur schema for root sequences allows CDS names that match `^[a-zA-Z0-9*_-]+$` ([schema-export-root-sequence.json#L25](https://github.com/nextstrain/augur/blob/0dc610704848c3a472f7715a1161277ef9b625a4/augur/data/schema-export-root-sequence.json#L25)), while the augur schema for annotations allows `^(?!nuc$)[a-zA-Z0-9*_.()-]+$`, which adds `.`, `(` and `)` ([schema-annotations.json#L31](https://github.com/nextstrain/augur/blob/0dc610704848c3a472f7715a1161277ef9b625a4/augur/data/schema-annotations.json#L31)). A CDS name such as `ORF1a.1` passes the annotation schema and fails the root-sequence schema. The two patterns disagree inside augur; TreeTime takes CDS names from the annotation without a check.

## Fix direction

Check each CDS name against the root-sequence pattern when the outputs include amino-acid sequences keyed by CDS name, and report each name that fails in the run record ([M-cli-run-warnings-and-diagnostics-only-on-stderr.md](M-cli-run-warnings-and-diagnostics-only-on-stderr.md)).

> [!IMPORTANT]
> **Investigation required.** Check which TreeTime outputs carry CDS names in places that the root-sequence schema validates, today and after TreeTime writes the root's `aa_sequences` and `reference` ([kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md), axis 12).
