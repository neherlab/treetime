import type { RunRecord, RunResults, AncestralResults, BranchMutations, RecurrentSite } from "@neherlab/app-contracts";
import { useMemo } from "react";
import ChartColumn from "~icons/lucide/chart-column";

import { DataTable, dataColumns } from "../components/DataTable";
import { Panel, runTimeEntry, SummaryStrip, type SummaryEntry } from "../components/Panel";
import { Button } from "../ui/button";
import { OutputFiles } from "./OutputFiles";
import { settingText } from "./settingText";
import { MissingTree, TreeView, type TreeData } from "./TreeView";
import { useRerun } from "./useRerun";

const branchColumn = dataColumns<BranchMutations>();

const BRANCH_COLUMNS = branchColumn.columns([
  branchColumn.accessor((row) => row.name, {
    id: "branch",
    header: "Branch above",
    cell: ({ row }) =>
      row.original.tips === 1 ? row.original.name : `${row.original.name} (${row.original.tips} samples)`,
  }),
  branchColumn.accessor((row) => row.mutations.length, { id: "count", header: "Mutations" }),
  branchColumn.accessor((row) => row.mutations.join(" "), {
    id: "list",
    header: "List",
    cell: ({ getValue }) => <span className="font-mono whitespace-normal">{getValue()}</span>,
  }),
]);

const siteColumn = dataColumns<RecurrentSite>();

const SITE_COLUMNS = siteColumn.columns([
  siteColumn.accessor((row) => row.position, { id: "position", header: "Position" }),
  siteColumn.accessor((row) => row.branches, { id: "branches", header: "Branches" }),
]);

const BRANCH_NUMERIC = new Set(["count"]);

const SITE_NUMERIC = new Set(["position", "branches"]);

const COUNT_SORT = [{ id: "count", desc: true }];

const BRANCHES_SORT = [{ id: "branches", desc: true }];

export function AncestralResultsView({
  record,
  results,
  data,
  tree,
}: {
  record: RunRecord;
  results: RunResults;
  data: AncestralResults;
  tree: TreeData | undefined;
}) {
  const branches = data.branches;
  const sites = data.recurrent_sites;
  const total = data.mutations;
  const analyzeHomoplasy = useRerun(record, "homoplasy");

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
    <div className="grid gap-4">
      <SummaryStrip entries={summary} />
      {tree === undefined ? <MissingTree /> : <TreeView data={tree} colorBy={undefined} />}
      <div className="grid gap-4 @4xl:grid-cols-[minmax(0,2fr)_minmax(0,1fr)]">
        <Panel title="Branches with the most mutations">
          <DataTable
            label="Branches with the most mutations"
            columns={BRANCH_COLUMNS}
            rows={branches}
            rowId={branchKey}
            initialSorting={COUNT_SORT}
            numeric={BRANCH_NUMERIC}
          />
        </Panel>
        <Panel
          title="Sites mutated on several branches"
          actions={
            sites.length === 0 ? undefined : (
              <Button type="button" variant="outline" size="xs" onClick={analyzeHomoplasy}>
                <ChartColumn aria-hidden />
                Analyze homoplasy
              </Button>
            )
          }
        >
          <DataTable
            label="Sites mutated on several branches"
            columns={SITE_COLUMNS}
            rows={sites}
            rowId={siteKey}
            initialSorting={BRANCHES_SORT}
            numeric={SITE_NUMERIC}
          />
        </Panel>
      </div>
      <OutputFiles record={record} citation={results.citation} />
    </div>
  );
}

function branchKey(row: BranchMutations): string {
  return row.name;
}

function siteKey(row: RecurrentSite): string {
  return String(row.position);
}
