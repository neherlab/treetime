import type { RunRecord, AncestralResults, BranchMutations, RecurrentSite } from "@neherlab/app-contracts";
import { CellText, DataTable, dataColumns } from "@neherlab/design/ui/data-table";
import { useCallback, useMemo, type ReactNode } from "react";

import { GENOTYPE_COLORING, genotypeColorBy } from "../auspice/toolbar";
import { Panel, runTimeEntry, SummaryStrip, type SummaryEntry } from "../components/Panel";
import { useLocale } from "../i18n/locale";
import { settingText } from "./settingText";
import { MissingTree, TreeView, type TreeData } from "./TreeView";
import type { TreeLink } from "./TreeWorkspace";

const branchColumn = dataColumns<BranchMutations>();

const siteColumn = dataColumns<RecurrentSite>();

const COUNT_SORT = [{ id: "count", desc: true }];

const BRANCHES_SORT = [{ id: "branches", desc: true }];

export function AncestralResultsView({
  record,
  data,
  tree,
}: {
  record: RunRecord;
  data: AncestralResults;
  tree: TreeData | undefined;
}) {
  const branches = data.branches;
  const sites = data.recurrent_sites;
  const total = data.mutations;
  const locale = useLocale();

  const summary = useMemo<SummaryEntry[]>(
    () => [
      { label: "Mutations", value: String(total), detail: "On all branches" },
      { label: "Branches with mutations", value: String(branches.length) },
      { label: "Homoplasic sites", value: String(sites.length), detail: "Sites mutated on several branches" },
      { label: "Substitution model", value: settingText(record, "model", locale) },
      { label: "Reconstruction", value: settingText(record, "method_anc", locale) },
      runTimeEntry(record),
    ],
    [branches.length, locale, record, sites.length, total],
  );

  const colorBy = tree === undefined ? undefined : mutationColoring(tree, sites);

  const below = useCallback(
    (link: TreeLink | undefined) => <MutationTables branches={branches} sites={sites} link={link} />,
    [branches, sites],
  );

  return (
    <div className="grid gap-4">
      <SummaryStrip entries={summary} />
      {tree === undefined ? (
        <>
          <MissingTree />
          {below(undefined)}
        </>
      ) : (
        <TreeView data={tree} colorBy={colorBy} entropy={false} below={below} />
      )}
    </div>
  );
}

function MutationTables({
  branches,
  sites,
  link,
}: {
  branches: readonly BranchMutations[];
  sites: readonly RecurrentSite[];
  link: TreeLink | undefined;
}) {
  const select = link?.select;
  const setColorBy = link?.setColorBy;

  const onSite = useMemo(
    () => (setColorBy === undefined ? undefined : (position: string) => setColorBy(genotypeColorBy(Number(position)))),
    [setColorBy],
  );

  const branchTable = useMemo(() => branchColumns(select), [select]);
  const siteTable = useMemo(() => siteColumns(onSite), [onSite]);

  return (
    <div className="grid gap-4 @4xl:grid-cols-[minmax(0,2fr)_minmax(0,1fr)]">
      <Panel
        title="Mutations by branch"
        hint={link === undefined ? undefined : "Click a branch to show it in the tree"}
      >
        <DataTable
          label="Mutations by branch"
          columns={branchTable}
          rows={branches}
          rowId={branchKey}
          initialSorting={COUNT_SORT}
        />
      </Panel>
      <Panel
        title="Sites mutated on several branches"
        hint={link === undefined ? undefined : "Click a site to color the tree by its states"}
      >
        <DataTable
          label="Sites mutated on several branches"
          columns={siteTable}
          rows={sites}
          rowId={siteKey}
          initialSorting={BRANCHES_SORT}
        />
      </Panel>
    </div>
  );
}

function branchColumns(onBranch: ((name: string) => void) | undefined) {
  return branchColumn.columns([
    branchColumn.accessor((row) => row.name, {
      id: "branch",
      header: "Branch above",
      meta: { width: "1fr", minWidth: 160 },
      cell: ({ row }) => {
        const text =
          row.original.tips === 1 ? row.original.name : `${row.original.name} (${row.original.tips} samples)`;

        return onBranch === undefined ? (
          <CellText text={text} />
        ) : (
          <RowButton value={row.original.name} onClick={onBranch} label={`Show ${row.original.name} in the tree`}>
            {text}
          </RowButton>
        );
      },
    }),
    branchColumn.accessor((row) => row.mutations.length, {
      id: "count",
      header: "Mutations",
      meta: { width: 100, numeric: true },
    }),
    branchColumn.accessor((row) => row.mutations.join(" "), {
      id: "list",
      header: "List",
      meta: { width: "2fr", minWidth: 200 },
      cell: ({ getValue }) => <span className="font-mono whitespace-normal">{getValue()}</span>,
    }),
  ]);
}

function siteColumns(onSite: ((position: string) => void) | undefined) {
  return siteColumn.columns([
    siteColumn.accessor((row) => row.position, {
      id: "position",
      header: "Position",
      meta: { width: "1fr", numeric: true },
      cell: ({ row }) =>
        onSite === undefined ? (
          row.original.position
        ) : (
          <RowButton
            value={String(row.original.position)}
            onClick={onSite}
            label={`Color the tree by position ${row.original.position}`}
          >
            {row.original.position}
          </RowButton>
        ),
    }),
    siteColumn.accessor((row) => row.branches, {
      id: "branches",
      header: "Branches",
      meta: { width: "1fr", numeric: true },
    }),
  ]);
}

function RowButton({
  value,
  onClick,
  label,
  children,
}: {
  value: string;
  onClick: (value: string) => void;
  label: string;
  children: ReactNode;
}) {
  const click = useCallback(() => onClick(value), [onClick, value]);

  return (
    <button
      type="button"
      onClick={click}
      title={label}
      className="text-primary min-w-0 truncate text-left underline-offset-4 hover:underline"
    >
      {children}
    </button>
  );
}

function mutationColoring(tree: TreeData, sites: readonly RecurrentSite[]): string | undefined {
  const top = sites[0];
  const genotype = tree.tree.colorings.some((coloring) => coloring.key === GENOTYPE_COLORING);

  return top === undefined || !genotype ? undefined : genotypeColorBy(top.position);
}

function branchKey(row: BranchMutations): string {
  return row.name;
}

function siteKey(row: RecurrentSite): string {
  return String(row.position);
}
