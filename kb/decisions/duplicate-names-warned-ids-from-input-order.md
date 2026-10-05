# Duplicate names are warned about, and identity comes from input order

Sample and node names repeat in real-world input. TreeTime treats such input as degraded: it warns where detection is cheap, runs as if names were unique, and never uses names for internal identity.

**Type**: Design rule for v1.

**Status**: Approved on 2026-10-06.

## Rule

- **Identity**: code identifies every entry by an ID derived from its input order: the order of records in a FASTA file, the order of nodes in a Newick tree. Code never identifies an entry by its name
- **Names**: names are metadata. They pass from input to output unchanged, through a mapping from ID to name, and have no internal use. Outputs and screen reports restore names from IDs where they need them
- **Detection**: where a format or an algorithm can detect duplicates with no or small cost in time and memory, it does, and every surface (CLI, web app, desktop app) shows a warning in its own form
- **Processing**: algorithms assume unique names, except where they can detect duplicates cheaply and avoid the corruption that duplicates cause
- **No guessing**: TreeTime does not guess which entry a duplicate name was meant to identify

## Why

- Duplicates are a data error of the user, but they are inevitable in real-world data, so TreeTime expects them and handles them where it can
- Data with duplicates is less reliable: duplicates often show a lack of attention to detail and come with other data problems, so the user must see a warning
- TreeTime can avoid only part of the damage. Identity by input order keeps a duplicate from corrupting the computation itself; formats keyed by name, such as augur node data JSON, cannot hold two entries with the same name

## v0 behavior

v0 has no duplicate check. `TreeAnc` looks leaves up in a dictionary keyed by name, which keeps the last leaf of each name (`_leaves_lookup`, [packages/legacy/treetime/treetime/treeanc.py#L459](../../packages/legacy/treetime/treetime/treeanc.py#L459)).

## Related

- [kb/issues/M-io-duplicate-node-names-overwrite-name-keyed-outputs.md](../issues/M-io-duplicate-node-names-overwrite-name-keyed-outputs.md)
- [kb/issues/M-io-sequence-name-matching-unreliable.md](../issues/M-io-sequence-name-matching-unreliable.md)
- [kb/issues/N-reroot-duplicated-tip-name-resolution.md](../issues/N-reroot-duplicated-tip-name-resolution.md)
- [kb/proposals/input-name-matching-validation.md](../proposals/input-name-matching-validation.md)
