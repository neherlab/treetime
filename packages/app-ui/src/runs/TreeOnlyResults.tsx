import type { RunRecordResult, RunResultsResult } from "@neherlab/app-contracts";
import { useMemo } from "react";

import { formatDuration } from "../format";
import type { TreeSummary } from "../results/types";
import { OutputFiles } from "./OutputFiles";
import { SummaryStrip, type SummaryEntry } from "./Panel";
import { MissingTree } from "./TimetreeResults";
import { TreeView, type TreeData } from "./TreeView";

export function TreeOnlyResults({
  record,
  results,
  data,
  tree,
}: {
  record: RunRecordResult;
  results: RunResultsResult;
  data: TreeSummary;
  tree: TreeData | undefined;
}) {
  const summary = useMemo<SummaryEntry[]>(() => {
    const length = data.total_branch_length ?? undefined;
    const model = data.substitution_model ?? undefined;

    return [
      { label: "Samples", value: tree === undefined ? "-" : String(data.samples) },
      { label: "Internal nodes", value: tree === undefined ? "-" : String(data.internal_nodes) },
      {
        label: "Mutations",
        value: tree === undefined ? "-" : String(data.mutations),
        detail: "On all branches of the written tree",
      },
      ...(length === undefined
        ? []
        : [
            {
              label: "Total branch length",
              value: length.toPrecision(4),
              detail: "Sum of the branch lengths in the node data",
            },
          ]),
      ...(model === undefined
        ? []
        : [
            {
              label: "Substitution model",
              value: model.name,
              detail: `Overall rate μ = ${model.mu.toPrecision(4)}`,
            },
          ]),
      {
        label: "Run time",
        value:
          record.duration_seconds === null || record.duration_seconds === undefined
            ? "-"
            : formatDuration(record.duration_seconds),
      },
    ];
  }, [data, record, tree]);

  return (
    <div className="grid gap-3.5">
      <SummaryStrip entries={summary} />
      {tree === undefined ? <MissingTree /> : <TreeView data={tree} colorBy={undefined} />}
      <OutputFiles record={record} methods={undefined} citation={results.citation} />
    </div>
  );
}
