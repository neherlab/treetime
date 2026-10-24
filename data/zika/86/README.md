# `zika/86`

86 Zika virus genomes sampled from 2013 to 2016, the Zika example of TreeTime 0.x.

## Metadata

`metadata.tsv` holds the TreeTime 0.x columns `#name`, `date`, and `country`, and a `region` column: the region of each country by the geolocation rules of augur (`augur/data/geolocation_rules.tsv`), with the Americas split into `north_america` and `south_america` as in Nextstrain metadata. Suriname and Cuba have no rule there and take the region of their neighbors.

The subsample `zika/20` takes its metadata rows from this table.

## Other files

- `aln.phylip.xz`: the sequences of `aln.fasta.xz` in the PHYLIP format of TreeTime 0.x
