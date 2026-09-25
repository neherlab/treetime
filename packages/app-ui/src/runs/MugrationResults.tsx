import type { RunRecordResult } from "@neherlab/app-contracts";
import { useMemo } from "react";

import { formatDuration } from "../format";
import type { RunResults } from "../results/load";
import {
  UNCERTAIN_STATE_PROBABILITY,
  ancestorStates,
  stateChanges,
  uncertainAncestors,
  type AncestorState,
  type StateChange,
} from "../results/mutations";
import { OutputFiles } from "./OutputFiles";
import { Panel, SummaryStrip, type SummaryEntry } from "./Panel";
import { settingText } from "./settingText";
import { SortableTable, type Column } from "./SortableTable";
import { MissingTree } from "./TimetreeResults";
import { TreeView } from "./TreeView";

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

export function MugrationResults({ record, results }: { record: RunRecordResult; results: RunResults }) {
  const attribute = results.traits?.attribute ?? settingText(record, "attribute");
  const tree = results.auspice?.tree;
  const changes = useMemo(() => (tree === undefined ? [] : stateChanges(tree, attribute)), [attribute, tree]);
  const uncertain = useMemo(() => (tree === undefined ? [] : uncertainAncestors(tree, attribute)), [attribute, tree]);

  const root = useMemo(
    () =>
      tree === undefined ? undefined : ancestorStates(tree, attribute).find((state) => state.name === tree.root.name),
    [attribute, tree],
  );

  const states = useMemo(
    () => new Set(tree?.tips.flatMap((tip) => tip.traits.get(attribute)?.value ?? []) ?? []),
    [attribute, tree],
  );

  const auspice = results.auspice;

  const summary = useMemo<SummaryEntry[]>(
    () => [
      { label: "Trait", value: attribute },
      { label: "States", value: String(states.size), detail: "Among the samples" },
      {
        label: "State changes",
        value: String(changes.reduce((sum, change) => sum + change.branches, 0)),
        detail: "Branches whose state differs from the parent",
      },
      {
        label: "Uncertain ancestors",
        value: String(uncertain.length),
        detail: `Most probable state below P = ${UNCERTAIN_STATE_PROBABILITY}`,
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
    [attribute, changes, record, root, states.size, uncertain.length],
  );

  return (
    <div className="grid gap-3.5">
      <SummaryStrip entries={summary} />
      {auspice === undefined ? <MissingTree /> : <TreeView auspice={auspice} colorBy={attribute} />}
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
        <Panel
          title="Uncertain ancestors"
          hint={`Ancestors whose most probable state has P below ${UNCERTAIN_STATE_PROBABILITY}`}
        >
          {uncertain.length === 0 ? (
            <p className="text-ink-muted px-3.5 py-3">
              Every ancestor has a state with P of at least {UNCERTAIN_STATE_PROBABILITY}.
            </p>
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
      <OutputFiles record={record} methods={undefined} />
    </div>
  );
}

function changeKey(row: StateChange): string {
  return JSON.stringify([row.from, row.to]);
}

function ancestorKey(row: AncestorState): string {
  return row.name;
}
