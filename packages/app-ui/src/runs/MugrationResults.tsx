import type { RunRecord, RunResults } from "@neherlab/app-contracts";
import { useMemo } from "react";

import { DataTable, dataColumns } from "../components/DataTable";
import { Panel, runTimeEntry, SummaryStrip, type SummaryEntry } from "../components/Panel";
import type { AncestorState, MugrationData, StateChange } from "../results/types";
import { Empty, EmptyDescription } from "../ui/empty";
import { OutputFiles } from "./OutputFiles";
import { MissingTree, TreeView, type TreeData } from "./TreeView";

const BRANCHES_SORT = [{ id: "branches", desc: true }];

const PROBABILITY_SORT = [{ id: "probability", desc: false }];

const CHANGE_NUMERIC = new Set(["branches"]);

const ANCESTOR_NUMERIC = new Set(["probability"]);

type StateColors = ReadonlyMap<string, string>;

const changeColumn = dataColumns<StateChange>();

const ancestorColumn = dataColumns<AncestorState>();

function changeColumns(colors: StateColors) {
  return changeColumn.columns([
    changeColumn.accessor((row) => row.from, {
      id: "from",
      header: "From",
      cell: ({ getValue }) => <StateLabel state={getValue()} colors={colors} />,
    }),
    changeColumn.accessor((row) => row.to, {
      id: "to",
      header: "To",
      cell: ({ getValue }) => <StateLabel state={getValue()} colors={colors} />,
    }),
    changeColumn.accessor((row) => row.branches, { id: "branches", header: "Branches" }),
  ]);
}

function ancestorColumns(colors: StateColors) {
  return ancestorColumn.columns([
    ancestorColumn.accessor((row) => row.name, {
      id: "name",
      header: "Ancestor",
      cell: ({ row }) => `${row.original.name} (${row.original.tips} samples)`,
    }),
    ancestorColumn.accessor((row) => row.state, {
      id: "state",
      header: "Most probable state",
      cell: ({ getValue }) => <StateLabel state={getValue()} colors={colors} />,
    }),
    ancestorColumn.accessor((row) => row.probability, {
      id: "probability",
      header: "Probability",
      cell: ({ getValue }) => getValue().toFixed(3),
    }),
  ]);
}

export function MugrationResults({
  record,
  results,
  data,
  tree,
}: {
  record: RunRecord;
  results: RunResults;
  data: MugrationData;
  tree: TreeData | undefined;
}) {
  const attribute = data.attribute;
  const changes = data.state_changes;
  const uncertain = data.uncertain_ancestors;
  const threshold = data.uncertain_below;
  const root = data.root ?? undefined;
  const colorings = tree?.tree.colorings;

  const [changeCols, ancestorCols] = useMemo(() => {
    const scale = colorings?.find((coloring) => coloring.key === attribute)?.scale ?? [];
    const colors: StateColors = new Map(scale.map((entry) => [entry.state, entry.color]));

    return [changeColumns(colors), ancestorColumns(colors)] as const;
  }, [attribute, colorings]);

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

  return (
    <div className="grid gap-4">
      <SummaryStrip entries={summary} />
      {tree === undefined ? <MissingTree /> : <TreeView data={tree} colorBy={attribute} />}
      <div className="grid gap-4 @4xl:grid-cols-2">
        <Panel title="State changes" hint="Counted over branches of the tree">
          <DataTable
            label="State changes"
            columns={changeCols}
            rows={changes}
            rowId={changeKey}
            numeric={CHANGE_NUMERIC}
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
              numeric={ANCESTOR_NUMERIC}
              initialSorting={PROBABILITY_SORT}
            />
          )}
        </Panel>
      </div>
      <OutputFiles record={record} citation={results.citation} />
    </div>
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
