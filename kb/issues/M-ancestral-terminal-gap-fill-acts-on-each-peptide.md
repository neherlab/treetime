# Terminal-gap fill acts on each peptide instead of the genome

`--gap-fill only-terminal` turns the leading and trailing gaps of each sequence into unknown characters. The amino-acid route applies it to each tip peptide (`fn read_cds_translations()` [packages/app-commands/src/commands/ancestral/run.rs#L382-L413](../../packages/app-commands/src/commands/ancestral/run.rs#L382-L413)), while the nucleotide route applies it to the whole genome. A deletion that touches the start or end of a CDS therefore becomes `X` in the amino-acid route, while it stays an interior gap of the genome in the nucleotide route, and the two routes disagree at that codon.

On `mpox/clade-ii/2000` this explains most of the 4,915 internal-node residues where one route writes `X` or `-` and the other a different character: the CDS with many such differences all have tips with CDS-terminal gap runs (OPG068 1,562 tip codons, OPG015_dup 1,175, OPG082 1,014), 11,413 tip codons in total. Method: section "Measured disagreement of the two routes" of [kb/reports/codon-aware-ancestral-reconstruction.md](../reports/codon-aware-ancestral-reconstruction.md#measured-disagreement-of-the-two-routes).

> [!IMPORTANT]
> **Investigation required.** RSV L has 136 differences of this class on `rsv/a/2000` but only 10 tip codons with CDS-terminal gaps, so the terminal-gap fill does not explain them. Trace their cause.

## Fix direction

When TreeTime translates its own reconstruction ([M-ancestral-aa-outputs-disagree-with-translated-nucleotide-outputs.md](M-ancestral-aa-outputs-disagree-with-translated-nucleotide-outputs.md)), the gap fill acts once on the genome, and the peptides inherit its result. `--translations` then remains only for protein-only input, where no genome exists and the fill of each peptide is the correct meaning. Until then, apply the fill to a peptide position only when its codon lies in a terminal gap run of the tip's genome.

## Validation

A tip with a deletion that covers the first codons of one CDS and nothing at the genome ends gives `-` at those codons in both routes.
