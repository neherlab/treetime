# Timetree divergence-tree output is unimplemented

The timetree command writes tree formats through `EdgeTimetree::nwk_weight()`, which returns `time_length`. Newick and Nexus therefore already contain time-based branch lengths. The command does not also emit the v0 divergence tree whose edge lengths are substitutions per site.

V0 writes both identities by serializing inferred time lengths and then mutation lengths. V1 retains the divergence value in `EdgeTimetree.base.branch_length`, but `write_timetree_tree_outputs()` accepts one graph and one output map, with no approved naming or CLI contract for a second tree identity.

## Decision required

- Define how one command exposes two tree identities through the existing output-selection and per-file flag model.
- Define which selected formats apply to each identity and how default paths avoid collisions.
- Preserve current timetree output as time-based; add divergence output from the stored substitution-per-site length.

Create an implementation ticket after the multi-identity output contract is approved.

## Related

- [kb/proposals/output-format-selection.md](../proposals/output-format-selection.md)
