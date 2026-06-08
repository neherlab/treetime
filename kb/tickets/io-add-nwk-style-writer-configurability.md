# Add --nwk-style flag for Newick output dialect configurability

Add a `--nwk-style` CLI flag to all tree-outputting commands controlling which annotation dialect the Newick writer uses. Applies to both `.nwk` and `.nexus` output (both embed the same Newick string).

## Values

- `plain` (default): no annotations, no comments. Round-trips through every parser including `bio::io::newick`
- `annotated`: BEAST-style `[&key=value,...]`. Compatible with BEAST, FigTree, MrBayes consensus, IQ-TREE, DendroPy, BioPython
- `nhx`: NHX `[&&NHX:key=value:...]`. Compatible with Forester, ETE, DendroPy

## Default rationale

Default `plain` because:

- The current v1 reader (`bio::io::newick`) rejects annotations (see [M-io-newick-output-incompatible-with-reader.md](../issues/M-io-newick-output-incompatible-with-reader.md))
- `.nexus` provides no additional annotation capability over `.nwk`
- v1 has never shipped; no backward compatibility concern

Once the custom parser (see [io-replace-bio-newick-parser-with-custom-dialect-aware-parser.md](io-replace-bio-newick-parser-with-custom-dialect-aware-parser.md)) ships, the default can be reconsidered.

## Implementation

- Add `NwkStyle` enum to `commands/shared/output.rs` (`OutputArgs`)
- Pass `NwkStyle` through to `nwk_write_with` and `nex_write_with`
- When `plain`: pass empty `CommentProviders` (suppress all annotations)
- When `annotated`: pass the command's `CommentProviders` with BEAST-style formatting (current behavior)
- When `nhx`: pass `CommentProviders` with NHX formatting (colon-separated, `&&NHX` prefix)
- Annotation content is command-determined (mutations for ancestral, mutations+dates for timetree, states for mugration). The flag controls syntax and inclusion, not content

## Affected commands

All commands using `write_graph_files_with_options`: optimize, ancestral, timetree, mugration, clock, prune.

## Orthogonal concerns

- eNewick and Rich Newick are structural (DAG topology, extra colon fields), not comment conventions. These belong on a separate `--output-tree-format` axis when network algorithms ship
- Per-file output paths (`--output-tree-nwk`, `--output-tree-nexus`) are orthogonal to `--nwk-style` (WHERE vs HOW)

## Related issues

- Source: [M-io-newick-output-incompatible-with-reader.md](../issues/M-io-newick-output-incompatible-with-reader.md)
- Source: [N-io-nexus-writer-no-annotation-advantage.md](../issues/N-io-nexus-writer-no-annotation-advantage.md)
- [kb/reports/newick-annotation-dialects.md](../reports/newick-annotation-dialects.md)
