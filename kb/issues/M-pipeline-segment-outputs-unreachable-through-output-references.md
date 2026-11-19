# Pipeline output references do not reach the per-segment outputs of a step on segment trees

A `pipeline` step passes an output of an earlier step with `{{ steps.<name>.outputs.<selection> }}` (`fn resolve_selection_path` [packages/app-cli/src/cli/pipeline/resolve.rs](../../packages/app-cli/src/cli/pipeline/resolve.rs)). A step that reads a graph of segment trees (`segment_trees`, `segment_alignments`) writes its tree outputs once per segment, numbered in the order of the lists (`optimize.1.nexus`, `optimize.2.nexus`, ...), and no reference reaches these files:

- **Per-segment outputs are rejected**: `{{ steps.optimize.outputs.nexus }}` resolves to the four Nexus files of a step with four segments, and the pipeline stops with "resolves to one file per alignment file (...); reference one of these paths instead". The reference syntax has no way to name one of these files, and no way to pass all of them to a list such as `segment_trees` of the next step
- **Whole-run outputs resolve to files that are never written**: `{{ steps.optimize.outputs.nwk }}` resolves to `optimize.nwk`. A step on segment trees does not write that file ("Not writing ...: these files hold one tree for the whole run", `fn write_whole_graph_outputs` [packages/app-commands/src/commands/shared/segment_outputs.rs](../../packages/app-commands/src/commands/shared/segment_outputs.rs)), but the pipeline accepts the reference. In `segment_trees: ["{{ steps.optimize.outputs.nwk }}"]` of a timetree step, the run stops later with "they list 1 trees and 4 alignments", which does not name the reference as the cause. The resolution reads the planned outputs (`fn paths_by_selection` [packages/app-output/src/output_plan.rs](../../packages/app-output/src/output_plan.rs)), which do not know that the step reads segment trees
- **Plain Newick per segment needs a path**: a segment step writes Newick per segment only with `{alignment}` in `output_tree_nwk`. In a pipeline file this relative path resolves from the folder of the file, not from the output folder of the step, and `{{ output_all }}` is not a template namespace (only `env`, `steps`, `vars`)

Reproduced on `data/simulated/reassortment/80` with a two-step pipeline (`optimize`, then `timetree` on its trees).

The examples chain segment steps through `{{ steps.<name>.output_all }}` and the numbered file names, which the pipeline does not check before the run: `"{{ steps.optimize.output_all }}/optimize.1.nexus"` (`data/toy/reassortment/pipeline-optimized-graph.yaml`, `data/simulated/reassortment/80/pipeline-optimized-graph.yaml`).

> [!IMPORTANT]
> **Decision required.** How a reference names per-segment outputs. Options:
>
> - A whole-value reference to a per-segment output expands to the list of its files, in segment order: `segment_trees: "{{ steps.optimize.outputs.nexus }}"`. This matches how a whole-value reference to a list in `vars` keeps its type
> - An index selects one file: `{{ steps.optimize.outputs.nexus[1] }}`
>
> Either way, the planned outputs of a step on segment trees must list the files that the step writes, so that a reference to a whole-run file that the step skips fails at resolution with a message that names the reference.
