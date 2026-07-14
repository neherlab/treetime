# Auspice v2 output omits required updated metadata

The shared TreeIR Auspice writer declares schema version `v2` but omits required `meta.updated`. Every command using the shared writer can emit schema-invalid JSON.

The pinned Augur export-v2 schema requires `updated` [[spec](https://github.com/nextstrain/augur/blob/d8e38736037ba9474a809f9a5a63bc2b279d2407/augur/data/schema-export-v2.json#L14-L22)], and Augur populates it during export [[src](https://github.com/nextstrain/augur/blob/d8e38736037ba9474a809f9a5a63bc2b279d2407/augur/export_v2.py#L1159-L1162)].

## Potential solutions

- O1. Restore a required typed `updated` field and validate every generated document against the pinned schema.
- O2. Downgrade or change the declared schema contract. This would diverge from the selected Auspice v2 output and requires explicit approval.

## Recommendation

Require the output plan to generate one UTC calendar date in `YYYY-MM-DD` form and pass it into the writer as typed metadata. Validate ancestral, mugration, and timetree documents against the pinned official schema. Tests pass a fixed generation date; production supplies the current UTC date once per command run.

## Related issues

- [M-timetree-tree-output-inference-metadata-incomplete.md](M-timetree-tree-output-inference-metadata-incomplete.md)
