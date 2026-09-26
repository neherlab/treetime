import type { RunRecord, RunResults } from "@neherlab/app-contracts";
import { useMemo } from "react";

import type { AncestralData, BranchMutations, RecurrentSite } from "../results/types";
import { OutputFiles } from "./OutputFiles";
import { Panel, runTimeEntry, SummaryStrip, type SummaryEntry } from "./Panel";
import { settingText } from "./settingText";
import { SortableTable, type Column } from "./SortableTable";
import { MissingTree, TreeView, type TreeData } from "./TreeView";

const COUNT_SORT = { key: "count", descending: true };

const BRANCHES_SORT = { key: "branches", descending: true };

const BRANCH_COLUMNS: ReadonlyArray<Column<BranchMutations>> = [
  {
    key: "branch",
    label: "Branch above",
    kind: "text",
    value: (row) => row.name,
    render: (row) => (row.tips === 1 ? row.name : `${row.name} (${row.tips} samples)`),
  },
  { key: "count", label: "Mutations", kind: "number", value: (row) => row.mutations.length },
  {
    key: "list",
    label: "List",
    kind: "text",
    value: (row) => row.mutations.join(" "),
    render: (row) => <span className="font-mono whitespace-normal">{row.mutations.join(" ")}</span>,
  },
];

const SITE_COLUMNS: ReadonlyArray<Column<RecurrentSite>> = [
  { key: "position", label: "Position", kind: "number", value: (row) => row.position },
  { key: "branches", label: "Branches", kind: "number", value: (row) => row.branches },
];

export function AncestralResults({
  record,
  results,
  data,
  tree,
}: {
  record: RunRecord;
  results: RunResults;
  data: AncestralData;
  tree: TreeData | undefined;
}) {
  const branches = data.branches;
  const sites = data.recurrent_sites;
  const total = data.mutations;

  const summary = useMemo<SummaryEntry[]>(
    () => [
      { label: "Mutations", value: String(total), detail: "On all branches" },
      { label: "Branches with mutations", value: String(branches.length) },
      { label: "Sites mutated on several branches", value: String(sites.length), detail: "Candidates for homoplasy" },
      { label: "Substitution model", value: settingText(record, "model") },
      { label: "Reconstruction", value: settingText(record, "method_anc") },
      runTimeEntry(record),
    ],
    [branches.length, record, sites.length, total],
  );

  return (
    <div className="grid gap-3.5">
      <SummaryStrip entries={summary} />
      {tree === undefined ? <MissingTree /> : <TreeView data={tree} colorBy={undefined} />}
      <div className="grid gap-3.5 xl:grid-cols-[minmax(0,2fr)_minmax(0,1fr)]">
        <Panel title="Branches with the most mutations">
          <SortableTable
            label="Branches with the most mutations"
            columns={BRANCH_COLUMNS}
            rows={branches}
            rowKey={branchKey}
            initialSort={COUNT_SORT}
          />
        </Panel>
        <Panel title="Sites mutated on several branches">
          <SortableTable
            label="Sites mutated on several branches"
            columns={SITE_COLUMNS}
            rows={sites}
            rowKey={siteKey}
            initialSort={BRANCHES_SORT}
          />
        </Panel>
      </div>
      <OutputFiles record={record} methods={undefined} citation={results.citation} />
    </div>
  );
}

function branchKey(row: BranchMutations): string {
  return row.name;
}

function siteKey(row: RecurrentSite): string {
  return String(row.position);
}
