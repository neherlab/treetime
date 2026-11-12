import type { RunRecord, TreeChange, TreeStatistics } from "@neherlab/app-contracts";
import { useMemo } from "react";

import { runTimeEntry, SummaryStrip, type SummaryEntry } from "../components/Panel";
import { MissingTree, TreeView, type TreeData } from "./TreeView";

export function TreeOnlyResults({
  record,
  data,
  tree,
}: {
  record: RunRecord;
  data: TreeChange;
  tree: TreeData | undefined;
}) {
  const summary = useMemo<SummaryEntry[]>(() => {
    const { before, after } = data;

    const change = (
      label: string,
      value: (statistics: TreeStatistics) => number | undefined,
      format: (value: number) => string,
    ): SummaryEntry => {
      const was = before === undefined ? undefined : value(before);
      const now = after === undefined ? undefined : value(after);
      const entry: SummaryEntry = { label, value: now === undefined ? "-" : format(now) };

      if (was !== undefined) {
        entry.detail = `Input tree: ${format(was)}`;
      }

      return entry;
    };

    return [
      change("Samples", (statistics) => statistics.samples, String),
      change("Internal nodes", (statistics) => statistics.internal_nodes, String),
      change("Polytomies", (statistics) => statistics.polytomies, String),
      {
        label: "Mutations",
        value: after === undefined ? "-" : String(data.mutations),
        detail: "On all branches of the written tree",
      },
      change(
        "Total branch length",
        (statistics) => statistics.total_branch_length ?? undefined,
        (value) => value.toPrecision(4),
      ),
      runTimeEntry(record),
    ];
  }, [data, record]);

  return (
    <div className="grid gap-4">
      <SummaryStrip entries={summary} />
      {tree === undefined ? <MissingTree /> : <TreeView data={tree} colorBy={undefined} />}
    </div>
  );
}
