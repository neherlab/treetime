# Timetree divergence-tree output is unimplemented

> [!IMPORTANT]
> **Decision required.** The output contract for two tree identities in one command (naming, CLI flags, format selection, and default paths) is not approved. Open points:
>
> - How one command exposes two tree identities through the existing output-selection and per-file flag model
> - Which selected formats apply to each identity, and how default paths avoid collisions
>
> The current timetree output stays time-based; the divergence output comes from the stored substitution-per-site lengths.

`fn write_timetree_trees()` [packages/app-commands/src/commands/timetree/run.rs#L342-L377](../../packages/app-commands/src/commands/timetree/run.rs#L342-L377) builds one `AnnotatedGraph` with both length sets: `divergence_branch_lengths` (substitutions per site) and `time_branch_lengths` (years). The Newick and Nexus writers read `AnnotatedGraph::tree_branch_lengths()` [packages/app-output/src/annotated_graph.rs#L51-L55](../../packages/app-output/src/annotated_graph.rs#L51-L55), which returns the time lengths when they are set. Newick and Nexus therefore contain time-based branch lengths. The divergence lengths reach only the DOT output, the Auspice `div` values, and the `mutation_length` of the augur node data. The command does not emit the v0 divergence tree in Newick or Nexus.

V0 writes both identities by serializing inferred time lengths and then mutation lengths. V1 has the divergence lengths at the writer, but one `AnnotatedGraph` maps to one file per tree format, with no approved naming or CLI contract for a second tree identity.

Start implementation after the multi-identity output contract is approved.

## Related

- [kb/proposals/output-format-selection.md](../proposals/output-format-selection.md)
