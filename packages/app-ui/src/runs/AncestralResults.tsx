import type { RunRecordResult } from "@neherlab/app-contracts";
import { useMemo } from "react";

import { formatDuration } from "../format";
import type { RunResults } from "../results/load";
import {
  branchesByMutationCount,
  recurrentSites,
  type BranchMutations,
  type RecurrentSite,
} from "../results/mutations";
import { OutputFiles } from "./OutputFiles";
import { Panel, SummaryStrip, type SummaryEntry } from "./Panel";
import { settingText } from "./settingText";
import { SortableTable, type Column } from "./SortableTable";
import { MissingTree } from "./TimetreeResults";
import { TreeView } from "./TreeView";

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

export function AncestralResults({ record, results }: { record: RunRecordResult; results: RunResults }) {
  const tree = results.auspice?.tree;
  const branches = useMemo(() => (tree === undefined ? [] : branchesByMutationCount(tree)), [tree]);
  const sites = useMemo(() => (tree === undefined ? [] : recurrentSites(tree)), [tree]);
  const total = branches.reduce((sum, branch) => sum + branch.mutations.length, 0);
  const auspice = results.auspice;

  const summary = useMemo<SummaryEntry[]>(
    () => [
      { label: "Mutations", value: String(total), detail: "On all branches" },
      { label: "Branches with mutations", value: String(branches.length) },
      { label: "Sites mutated on several branches", value: String(sites.length), detail: "Candidates for homoplasy" },
      { label: "Substitution model", value: settingText(record, "model") },
      { label: "Reconstruction", value: settingText(record, "method_anc") },
      {
        label: "Run time",
        value:
          record.duration_seconds === null || record.duration_seconds === undefined
            ? "-"
            : formatDuration(record.duration_seconds),
      },
    ],
    [branches.length, record, sites.length, total],
  );

  return (
    <div className="grid gap-3.5">
      <SummaryStrip entries={summary} />
      {auspice === undefined ? <MissingTree /> : <TreeView auspice={auspice} colorBy={undefined} />}
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
      <OutputFiles record={record} methods={undefined} />
    </div>
  );
}

function branchKey(row: BranchMutations): string {
  return row.name;
}

function siteKey(row: RecurrentSite): string {
  return String(row.position);
}
