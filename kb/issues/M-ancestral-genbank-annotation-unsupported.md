# GenBank annotation format not supported for AA reconstruction

The `--annotation` flag in `treetime ancestral` accepts GFF3 format only. GenBank flat-file format (`.gb`, `.gbk`) is not parsed.

## Impact

The ncov (SARS-CoV-2) nextstrain build uses a GenBank annotation file (`config/reference_seq.gb`). Users of ncov and other builds that rely on GenBank annotations must convert to GFF3 before running `treetime ancestral` with per-CDS AA reconstruction.

augur's `reconstruct_translations()` delegates to BioPython's `SeqIO.read(format="genbank")` for GenBank parsing. nextclade accepts GFF3 only.

## Workaround

Convert GenBank to GFF3 before invoking treetime. Tools: `bp_genbank2gff3` (BioPerl), `genometools gt gff3` after `gt convertseq`, or manual extraction of CDS features.

## Resolution path

Add a GenBank parser to `read_gff3_annotations()` (or a parallel `read_genbank_annotations()`) that extracts CDS features with coordinates, strand, and gene name. Detect format by file extension or content sniffing.
