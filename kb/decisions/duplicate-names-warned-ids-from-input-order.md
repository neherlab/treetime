# Duplicate names are warned about, and identity comes from input order

Sample and node names repeat in real-world input. TreeTime treats such input as degraded: it warns where detection is cheap, runs as if names were unique, and never uses names for internal identity.

**Type**: Design rule for v1.

**Status**: Approved on 2026-10-06.

## Rule

- **Identity**: code identifies every entry by an ID derived from its input order: the order of records in a FASTA file, the order of nodes in a Newick tree. Code never identifies an entry by its name
- **Names**: names are metadata. They pass from input to output unchanged, through a mapping from ID to name, and have no internal use. Outputs and screen reports restore names from IDs where they need them
- **Pairing across inputs**: inputs that can only be paired by name (tree leaves, FASTA records, metadata rows) pair with the first entry in input order when a name has more than one entry. Later entries with that name get nothing from the pairing
- **Names in arguments**: a name that the user types, for example a tip of `--reroot-tips`, resolves to the first matching node in input order
- **Detection cost**: TreeTime adds no search for duplicates that costs time or memory on every run. It warns only where a format, a pairing or an algorithm finds a duplicate at no or small cost, for example while it builds an index that it needs anyway
- **Warning form**: the run result carries each warning as a typed entry, so the web and desktop apps show it as a notice next to the results, and the run log and the CLI show the same text
- **Processing**: algorithms assume unique names, except where they can detect duplicates cheaply and avoid the corruption that duplicates cause
- **No guessing**: TreeTime does not try to find out which entry a duplicate name was meant to identify

## Why

- Duplicates are a data error of the user, but they are inevitable in real-world data, so TreeTime expects them and handles them where it can
- Data with duplicates is less reliable: duplicates often show a lack of attention to detail and come with other data problems, so the user must see a warning when TreeTime finds them
- Good data must not pay for bad data. An exhaustive search for duplicates costs every user time and memory, also users whose data has none. A user who gives data whose samples cannot be resolved uniquely bears the result: output that can be wrong or corrupted
- Identity by input order keeps a duplicate from corrupting the computation itself. Formats keyed by name, such as augur node data JSON, cannot hold two entries with the same name, so there the warning is the only remedy

## v0 behavior

v0 has no duplicate check. `TreeAnc` looks leaves up in a dictionary keyed by name, which keeps the last leaf of each name (`_leaves_lookup`, [packages/legacy/treetime/treetime/treeanc.py#L459](../../packages/legacy/treetime/treetime/treeanc.py#L459)).

## Related

- [kb/issues/M-io-duplicate-node-names-overwrite-name-keyed-outputs.md](../issues/M-io-duplicate-node-names-overwrite-name-keyed-outputs.md)
- [kb/issues/M-io-sequence-name-matching-unreliable.md](../issues/M-io-sequence-name-matching-unreliable.md)
- [kb/issues/N-reroot-duplicated-tip-name-resolution.md](../issues/N-reroot-duplicated-tip-name-resolution.md)
- [kb/proposals/input-name-matching-validation.md](../proposals/input-name-matching-validation.md)
