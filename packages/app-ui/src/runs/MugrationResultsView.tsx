import type { RunRecord, AncestorState, MugrationResults, StateChange } from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { CellText, DataTable, dataColumns } from "@neherlab/design/ui/data-table";
import { Empty, EmptyDescription } from "@neherlab/design/ui/empty";
import { useCallback, useMemo } from "react";

import { Panel, runTimeEntry, SummaryStrip, type SummaryEntry } from "../components/Panel";
import { MissingTree, TreeView, type TreeData } from "./TreeView";
import type { TreeLink } from "./TreeWorkspace";

const BRANCHES_SORT = [{ id: "branches", desc: true }];

const PROBABILITY_SORT = [{ id: "probability", desc: false }];

type StateColors = ReadonlyMap<string, string>;

const NO_COLORS: StateColors = new Map();

const changeColumn = dataColumns<StateChange>();

const ancestorColumn = dataColumns<AncestorState>();

export function MugrationResultsView({
  record,
  data,
  tree,
}: {
  record: RunRecord;
  data: MugrationResults;
  tree: TreeData | undefined;
}) {
  const attribute = data.attribute;
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
      runTimeEntry(record),
    ],
    [attribute, data, record, root, threshold, uncertain.length],
  );

  const below = useCallback((link: TreeLink) => <LinkedTables data={data} link={link} />, [data]);

  return (
    <div className="grid gap-4">
      <SummaryStrip entries={summary} />
      {tree === undefined ? (
        <>
          <MissingTree />
          <MugrationTables data={data} colors={NO_COLORS} onSelect={undefined} />
        </>
      ) : (
        <TreeView data={tree} colorBy={attribute} below={below} />
      )}
    </div>
  );
}

function LinkedTables({ data, link }: { data: MugrationResults; link: TreeLink }) {
  const { stateColors, select } = link;
  const colors = useMemo(() => stateColors(data.attribute), [data.attribute, stateColors]);

  return <MugrationTables data={data} colors={colors} onSelect={select} />;
}

function MugrationTables({
  data,
  colors,
  onSelect,
}: {
  data: MugrationResults;
  colors: StateColors;
  onSelect: ((name: string) => void) | undefined;
}) {
  const uncertain = data.uncertain_ancestors;
  const threshold = data.uncertain_below;

  const [changeCols, ancestorCols] = useMemo(
    () => [changeColumns(colors), ancestorColumns(colors, onSelect)] as const,
    [colors, onSelect],
  );

  return (
    <div className="grid gap-4 @4xl:grid-cols-2">
      <Panel title="State changes" hint="Counted over branches of the tree">
        <DataTable
          label="State changes"
          columns={changeCols}
          rows={data.state_changes}
          rowId={changeKey}
          initialSorting={BRANCHES_SORT}
        />
      </Panel>
      <Panel title="Uncertain ancestors" hint={`Ancestors whose most probable state has P below ${threshold}`}>
        {uncertain.length === 0 ? (
          <Empty className="py-6">
            <EmptyDescription>Every ancestor has a state with P of at least {threshold}.</EmptyDescription>
          </Empty>
        ) : (
          <DataTable
            label="Uncertain ancestors"
            columns={ancestorCols}
            rows={uncertain}
            rowId={ancestorKey}
            initialSorting={PROBABILITY_SORT}
          />
        )}
      </Panel>
    </div>
  );
}

function changeColumns(colors: StateColors) {
  return changeColumn.columns([
    changeColumn.accessor((row) => row.from, {
      id: "from",
      header: "From",
      meta: { width: "1fr", minWidth: 120 },
      cell: ({ getValue }) => <StateLabel state={getValue()} colors={colors} />,
    }),
    changeColumn.accessor((row) => row.to, {
      id: "to",
      header: "To",
      meta: { width: "1fr", minWidth: 120 },
      cell: ({ getValue }) => <StateLabel state={getValue()} colors={colors} />,
    }),
    changeColumn.accessor((row) => row.branches, {
      id: "branches",
      header: "Branches",
      meta: { width: 90, numeric: true },
    }),
  ]);
}

function ancestorColumns(colors: StateColors, onSelect: ((name: string) => void) | undefined) {
  return ancestorColumn.columns([
    ancestorColumn.accessor((row) => row.name, {
      id: "name",
      header: "Ancestor",
      meta: { width: "2fr", minWidth: 160 },
      cell: ({ row }) => <AncestorName ancestor={row.original} onSelect={onSelect} />,
    }),
    ancestorColumn.accessor((row) => row.state, {
      id: "state",
      header: "Most probable state",
      meta: { width: "1fr", minWidth: 140 },
      cell: ({ getValue }) => <StateLabel state={getValue()} colors={colors} />,
    }),
    ancestorColumn.accessor((row) => row.probability, {
      id: "probability",
      header: "Probability",
      meta: { width: 110, numeric: true },
      cell: ({ getValue }) => getValue().toFixed(3),
    }),
  ]);
}

function AncestorName({
  ancestor,
  onSelect,
}: {
  ancestor: AncestorState;
  onSelect: ((name: string) => void) | undefined;
}) {
  const text = `${ancestor.name} (${ancestor.tips} samples)`;
  const select = useCallback(() => onSelect?.(ancestor.name), [ancestor.name, onSelect]);

  if (onSelect === undefined) {
    return <CellText text={text} />;
  }

  return (
    <Button type="button" variant="ghost" size="xs" onClick={select}>
      {text}
    </Button>
  );
}

function StateLabel({ state, colors }: { state: string; colors: StateColors }) {
  const color = colors.get(state);

  return (
    <span className="inline-flex items-center gap-1.5">
      {color !== undefined && (
        <svg aria-hidden viewBox="0 0 10 10" className="size-2.5 shrink-0">
          <circle cx="5" cy="5" r="5" fill={color} />
        </svg>
      )}
      {state}
    </span>
  );
}

function changeKey(row: StateChange): string {
  return JSON.stringify([row.from, row.to]);
}

function ancestorKey(row: AncestorState): string {
  return row.name;
}
