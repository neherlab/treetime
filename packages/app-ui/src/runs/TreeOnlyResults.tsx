import type { RunRecord, TreeSummary } from "@neherlab/app-contracts";
import { useMemo } from "react";

import { runTimeEntry, SummaryStrip, type SummaryEntry } from "../components/Panel";
import { MissingTree, TreeView, type TreeData } from "./TreeView";

export function TreeOnlyResults({
  record,
  data,
  tree,
}: {
  record: RunRecord;
  data: TreeSummary;
  tree: TreeData | undefined;
}) {
  const summary = useMemo<SummaryEntry[]>(() => {
    const length = data.total_branch_length ?? undefined;

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
      runTimeEntry(record),
    ];
  }, [data, record, tree]);

  return (
    <div className="grid gap-4">
      <SummaryStrip entries={summary} />
      {tree === undefined ? <MissingTree /> : <TreeView data={tree} colorBy={undefined} />}
    </div>
  );
}
