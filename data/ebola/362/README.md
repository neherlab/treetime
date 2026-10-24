# `ebola/362`

Ebola virus genomes of the 2014-2016 West African epidemic, the Ebola example of TreeTime 0.x.

## Metadata

`metadata.tsv` holds the TreeTime 0.x columns `name` and `" date"`, and three columns taken from the descriptions of `aln.fasta.xz`, which have the form `<name>|ebola|<accession>|<date>|<region>|<country>|<division>|<location>|...`:

- `region`: `africa`
- `country`: `guinea`, `liberia`, `sierra_leone`
- `division`: the district, 35 values

The subsamples `ebola/20` and `ebola/100` take their metadata rows from this table.
