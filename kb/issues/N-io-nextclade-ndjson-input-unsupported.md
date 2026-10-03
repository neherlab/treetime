# Nextclade NDJSON is not accepted as alignment input

Nextclade writes one result per sample with substitutions, deletion ranges and missing ranges relative to the dataset reference [[src](https://github.com/nextstrain/nextclade/blob/629ea8fb0cfe3165db13ddc75dab16d3c2292dbd/packages/nextclade/src/types/outputs.rs#L73-L86)]. Sparse partitions could take these lists directly, so a user who has run Nextclade would not need to write and read a full alignment.

No tool is known to read Nextclade NDJSON as phylogenetic input, so this use is untested.

## Open questions

- Is this input worth supporting next to MAPLE and VCF, which carry the same information in established formats?
- Nextclade output contains no reference sequence. Does the user pass the reference separately, or does TreeTime read it from the Nextclade dataset?
- How are Nextclade insertions handled, given the decision to treat gaps as missing data in MAT output?

## Related

- [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md): ranking of format work
- [N-io-maple-alignment-input-unsupported.md](N-io-maple-alignment-input-unsupported.md): another sparse alignment input
