# Knowledge base contains internal links with missing targets

A repository-wide Markdown target check finds broken internal links in the following tracked documents.

## Algorithm and test inventories

- `kb/algo/timetree.md`: three links to the removed `packages/treetime/src/commands/timetree/refinement.rs` path.
- `kb/tests/timetree.md`: links to removed `commands/timetree/output/auspice.rs` source and test paths.
- `kb/tests/mugration.md`: link to removed `packages/treetime/src/mugration/input.rs`.
- `kb/tests/supporting.md`: ten link occurrences covering seven distinct removed `packages/treetime-cli/src/convert/` source and test targets.

## Decisions and proposals

- `kb/decisions/command-optimize-standalone.md`: link to missing `N-timetree-node-data-confidence-not-emitted.md`.
- `kb/decisions/multi-format-tree-io.md`: five links to removed conversion and PhyloXML source paths.
- `kb/proposals/mugration-full-reconstruction-per-iteration.md`: link to removed `ancestral/gtr_inference_dense.rs`.
- `kb/proposals/output-format-selection.md`: five links to missing issues and tickets.

## Reports

- `kb/reports/augur-node-data-json.md`: four links to missing issues and tickets.
- `kb/reports/dense-openblas-profiling.md`: links to untracked `.build/docker/profiling/treetime` and `.out/treetime` artifacts.
- `kb/reports/iterative-tree-refinement/9-iteration-loop.md`: link to removed `commands/timetree/refinement.rs`.
- `kb/reports/ancestral-mugration-comparison/README.md`: link to missing `ancestral-iterative-gtr-refinement.md`.

The missing targets prevent readers and automated checks from following the KB's evidence graph. Each link must be updated to the current tracked artifact or removed when no equivalent exists.

## Related tickets

- [kb/tickets/kb-repair-missing-internal-link-targets.md](../tickets/kb-repair-missing-internal-link-targets.md)
