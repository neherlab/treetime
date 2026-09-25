import type { RunRecordResult } from "@neherlab/app-contracts";
import { useMemo } from "react";

import { formatDuration } from "../format";
import type { RunResults } from "../results/load";
import { OutputFiles } from "./OutputFiles";
import { SummaryStrip, type SummaryEntry } from "./Panel";
import { MissingTree } from "./TimetreeResults";
import { TreeView } from "./TreeView";

export function TreeOnlyResults({ record, results }: { record: RunRecordResult; results: RunResults }) {
  const tree = results.auspice?.tree;
  const mutations = tree?.nodes.reduce((sum, node) => sum + node.mutations.length, 0);

  const summary = useMemo<SummaryEntry[]>(
    () => [
      { label: "Samples", value: tree === undefined ? "-" : String(tree.tips.length) },
      { label: "Internal nodes", value: tree === undefined ? "-" : String(tree.nodes.length - tree.tips.length) },
      {
        label: "Mutations",
        value: mutations === undefined ? "-" : String(mutations),
        detail: "On all branches of the written tree",
      },
      ...(results.totalBranchLength === undefined
        ? []
        : [
            {
              label: "Total branch length",
              value: results.totalBranchLength.toPrecision(4),
              detail: "Sum of the branch lengths in the node data",
            },
          ]),
      ...(results.gtr === undefined
        ? []
        : [
            {
              label: "Substitution model",
              value: results.gtr.model,
              detail: `Overall rate μ = ${results.gtr.mu.toPrecision(4)}`,
            },
          ]),
      {
        label: "Run time",
        value:
          record.duration_seconds === null || record.duration_seconds === undefined
            ? "-"
            : formatDuration(record.duration_seconds),
      },
    ],
    [mutations, record, results.gtr, results.totalBranchLength, tree],
  );

  return (
    <div className="grid gap-3.5">
      <SummaryStrip entries={summary} />
      {results.auspice === undefined ? <MissingTree /> : <TreeView auspice={results.auspice} colorBy={undefined} />}
      <OutputFiles record={record} methods={undefined} />
    </div>
  );
}
