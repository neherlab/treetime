import type { RunRecord, RunResults } from "@neherlab/app-contracts";
import { useMemo } from "react";

import type { AncestorState, MugrationData, StateChange } from "../results/types";
import { OutputFiles } from "./OutputFiles";
import { Panel, runTimeEntry, SummaryStrip, type SummaryEntry } from "./Panel";
import { SortableTable, type Column } from "./SortableTable";
import { MissingTree, TreeView, type TreeData } from "./TreeView";

const BRANCHES_SORT = { key: "branches", descending: true };

const PROBABILITY_SORT = { key: "probability", descending: false };

type StateColors = ReadonlyMap<string, string>;

function changeColumns(colors: StateColors): ReadonlyArray<Column<StateChange>> {
  return [
    {
      key: "from",
      label: "From",
      kind: "text",
      value: (row) => row.from,
      render: (row) => <StateLabel state={row.from} colors={colors} />,
    },
    {
      key: "to",
      label: "To",
      kind: "text",
      value: (row) => row.to,
      render: (row) => <StateLabel state={row.to} colors={colors} />,
    },
    { key: "branches", label: "Branches", kind: "number", value: (row) => row.branches },
  ];
}

function ancestorColumns(colors: StateColors): ReadonlyArray<Column<AncestorState>> {
  return [
    {
      key: "name",
      label: "Ancestor",
      kind: "text",
      value: (row) => row.name,
      render: (row) => `${row.name} (${row.tips} samples)`,
    },
    {
      key: "state",
      label: "Most probable state",
      kind: "text",
      value: (row) => row.state,
      render: (row) => <StateLabel state={row.state} colors={colors} />,
    },
    {
      key: "probability",
      label: "Probability",
      kind: "number",
      value: (row) => row.probability,
      render: (row) => row.probability.toFixed(3),
    },
  ];
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
    <div className="grid gap-3.5">
      <SummaryStrip entries={summary} />
      {tree === undefined ? <MissingTree /> : <TreeView data={tree} colorBy={attribute} />}
      <div className="grid gap-3.5 @4xl:grid-cols-2">
        <Panel title="State changes" hint="Counted over branches of the tree">
          <SortableTable
            label="State changes"
            columns={changeCols}
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
              columns={ancestorCols}
              rows={uncertain}
              rowKey={ancestorKey}
              initialSort={PROBABILITY_SORT}
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
