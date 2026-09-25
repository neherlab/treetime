import type { RunRecordResult, RunResultsResult } from "@neherlab/app-contracts";
import { useMemo } from "react";

import { formatDuration } from "../format";
import type { AncestorState, MugrationData, StateChange } from "../results/types";
import { OutputFiles } from "./OutputFiles";
import { Panel, SummaryStrip, type SummaryEntry } from "./Panel";
import { SortableTable, type Column } from "./SortableTable";
import { MissingTree } from "./TimetreeResults";
import { TreeView, type TreeData } from "./TreeView";

const BRANCHES_SORT = { key: "branches", descending: true };

const PROBABILITY_SORT = { key: "probability", descending: false };

const CHANGE_COLUMNS: ReadonlyArray<Column<StateChange>> = [
  { key: "from", label: "From", kind: "text", value: (row) => row.from },
  { key: "to", label: "To", kind: "text", value: (row) => row.to },
  { key: "branches", label: "Branches", kind: "number", value: (row) => row.branches },
];

const ANCESTOR_COLUMNS: ReadonlyArray<Column<AncestorState>> = [
  {
    key: "name",
    label: "Ancestor",
    kind: "text",
    value: (row) => row.name,
    render: (row) => `${row.name} (${row.tips} samples)`,
  },
  { key: "state", label: "Most probable state", kind: "text", value: (row) => row.state },
  {
    key: "probability",
    label: "Probability",
    kind: "number",
    value: (row) => row.probability,
    render: (row) => row.probability.toFixed(3),
  },
];

export function MugrationResults({
  record,
  results,
  data,
  tree,
}: {
  record: RunRecordResult;
  results: RunResultsResult;
  data: MugrationData;
  tree: TreeData | undefined;
}) {
  const attribute = data.attribute;
  const changes = data.state_changes;
  const uncertain = data.uncertain_ancestors;
  const threshold = data.uncertain_below;
  const root = data.root ?? undefined;

  const summary = useMemo<SummaryEntry[]>(
    () => [
      { label: "Trait", value: attribute },
      { label: "States", value: String(data.states), detail: "Among the samples" },
      {
        label: "State changes",
        value: String(data.changed_branches),
        detail: "Branches whose state differs from the parent",
      },
      {
        label: "Uncertain ancestors",
        value: String(uncertain.length),
        detail: `Most probable state below P = ${threshold}`,
        tone: uncertain.length === 0 ? undefined : "caution",
      },
      {
        label: "Root state",
        value: root?.state ?? "-",
        detail: root === undefined ? undefined : `P = ${root.probability.toFixed(3)}`,
      },
      {
        label: "Run time",
        value:
          record.duration_seconds === null || record.duration_seconds === undefined
            ? "-"
            : formatDuration(record.duration_seconds),
      },
    ],
    [attribute, data, record, root, threshold, uncertain.length],
  );

  return (
    <div className="grid gap-3.5">
      <SummaryStrip entries={summary} />
      {tree === undefined ? <MissingTree /> : <TreeView data={tree} colorBy={attribute} />}
      <div className="grid gap-3.5 xl:grid-cols-2">
        <Panel title="State changes" hint="Counted over branches of the tree">
          <SortableTable
            label="State changes"
            columns={CHANGE_COLUMNS}
            rows={changes}
            rowKey={changeKey}
            initialSort={BRANCHES_SORT}
          />
        </Panel>
        <Panel title="Uncertain ancestors" hint={`Ancestors whose most probable state has P below ${threshold}`}>
          {uncertain.length === 0 ? (
            <p className="text-ink-muted px-3.5 py-3">Every ancestor has a state with P of at least {threshold}.</p>
          ) : (
            <SortableTable
              label="Uncertain ancestors"
              columns={ANCESTOR_COLUMNS}
              rows={uncertain}
              rowKey={ancestorKey}
              initialSort={PROBABILITY_SORT}
            />
          )}
        </Panel>
      </div>
      <OutputFiles record={record} methods={undefined} citation={results.citation} />
    </div>
  );
}

function changeKey(row: StateChange): string {
  return JSON.stringify([row.from, row.to]);
}

function ancestorKey(row: AncestorState): string {
  return row.name;
}
