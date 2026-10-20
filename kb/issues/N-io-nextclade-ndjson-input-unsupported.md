# Nextclade NDJSON is not accepted as alignment input

Nextclade writes one result per sample with substitutions, deletion ranges and missing ranges relative to the dataset reference [[src](https://github.com/nextstrain/nextclade/blob/629ea8fb0cfe3165db13ddc75dab16d3c2292dbd/packages/nextclade/src/types/outputs.rs#L73-L86)]. Sparse partitions could take these lists directly, so a user who has run Nextclade would not need to write and read a full alignment.

The alignment input of TreeTime reads FASTA only: `fn read_sequence_source()` passes every `--alignment` path to `fasta_read_file()` [packages/app-commands/src/commands/shared/sequence_inputs.rs#L84-L92](../../packages/app-commands/src/commands/shared/sequence_inputs.rs#L84-L92).

No tool is known to read Nextclade NDJSON as phylogenetic input, so this use is untested.

## Open questions

> [!IMPORTANT]
> **Decision required.** It is not decided whether Nextclade NDJSON is in scope. MAPLE and VCF carry the same per-sample differences from a reference in established formats, and no other phylogenetic tool reads Nextclade NDJSON. The options are:
>
> - Support Nextclade NDJSON as a sparse alignment input
> - Support only MAPLE or VCF for this use and close this issue
>
> If the format is in scope, two design questions remain:
>
> - Reference: Nextclade output contains no reference sequence. The user passes the reference separately, or TreeTime reads it from the Nextclade dataset
> - Insertions: Nextclade reports insertions relative to the reference. TreeTime drops them, or represents them in a way that agrees with the treatment of gaps as missing data in MAT output

## Related

- [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md): ranking of format work
- [N-io-maple-alignment-input-unsupported.md](N-io-maple-alignment-input-unsupported.md): another sparse alignment input
