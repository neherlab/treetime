# Restore and schema-validate Auspice updated metadata

Make every shared-writer Auspice v2 document satisfy the official pinned schema.

## Required changes

- Add required `meta.updated` to the typed TreeIR Auspice model.
- Generate the current UTC calendar date once per command output plan, format it as `YYYY-MM-DD`, and pass it explicitly into the writer.
- Make the generation date injectable so tests use a fixed date without normalizing output after serialization.
- Validate complete ancestral, mugration, and timetree JSON against the pinned Augur export-v2 schema.

## Validation

- Official-schema success for every command using the writer.
- Negative fixture proving omission of `updated` fails validation.
- Fixed injected generation date for golden tests and a production-boundary test for UTC `YYYY-MM-DD` formatting.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/H-io-auspice-v2-required-updated-missing.md](../issues/H-io-auspice-v2-required-updated-missing.md)
