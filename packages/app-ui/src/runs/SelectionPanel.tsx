import { useNavigate } from "@tanstack/react-router";
import { useCallback, useMemo } from "react";

import { formatDecimalDate } from "../format";
import { useRunList, useRunTrees } from "../queries";
import { indexClades, type CladeIndex } from "../results/clades";
import { intervalWidthDays } from "../results/estimates";
import { tipCount } from "../results/mutations";
import type { ResultTree } from "../results/tree";
import { Button } from "../ui";
import { DateIntervals, type DateRow } from "./DateIntervals";
import { Panel } from "./Panel";
import type { TreeLink } from "./TreeWorkspace";

export function SelectionPanel({
  runId,
  title,
  tree,
  link,
}: {
  runId: string;
  title: string;
  tree: ResultTree;
  link: TreeLink;
}) {
  const node = link.focus ?? tree.root;
  const others = useOtherTimetrees(runId);
  const navigate = useNavigate();
  const clades = useMemo(() => indexClades(tree), [tree]);
  const key = clades.keyOf.get(node) ?? "";
  const open = useCallback((id: string) => void navigate({ to: "/runs/$id/results", params: { id } }), [navigate]);

  const rows = useMemo(() => {
    const current: DateRow[] =
      node.date === undefined
        ? []
        : [{ id: runId, label: title, date: node.date, interval: node.dateInterval, current: true }];

    return [
      ...current,
      ...others.flatMap(({ id, label, clades: other }) => {
        const match = other.nodeOf.get(key);

        return match?.date === undefined
          ? []
          : [{ id, label, date: match.date, interval: match.dateInterval, current: false }];
      }),
    ];
  }, [key, node, others, runId, title]);

  const missing = others.length - (rows.length - (node.date === undefined ? 0 : 1));
  const width = intervalWidthDays(node.dateInterval);
  const isTip = node.children.length === 0;

  return (
    <Panel
      title={isTip ? node.name : `Ancestor of ${tipCount(node)} samples`}
      hint={
        link.focus === undefined || link.focus === tree.root
          ? "Click a branch in the tree to zoom into its clade"
          : undefined
      }
      actions={
        link.zoomed ? (
          <Button type="button" variant="ghost" size="sm" onClick={link.reset}>
            Show the whole tree
          </Button>
        ) : undefined
      }
    >
      <dl className="m-0 grid grid-cols-[auto_1fr] gap-x-3 gap-y-1 px-3.5 py-3 text-sm">
        <dt className="text-ink-faint">{isTip ? "Date in tree" : "Date"}</dt>
        <dd className="m-0 tabular-nums">{node.date === undefined ? "not dated" : formatDecimalDate(node.date)}</dd>
        {node.dateInterval !== undefined && width !== undefined && width > 0 && (
          <>
            <dt className="text-ink-faint">90% interval</dt>
            <dd className="m-0 tabular-nums">
              {formatDecimalDate(node.dateInterval[0])} to {formatDecimalDate(node.dateInterval[1])} (
              {Math.round(width)} days)
            </dd>
          </>
        )}
        {node.excluded === true && isTip && (
          <>
            <dt className="text-ink-faint">Clock</dt>
            <dd className="text-signal-warn m-0 font-bold">Excluded: no usable date or clock outlier</dd>
          </>
        )}
        <dt className="text-ink-faint">Mutations</dt>
        <dd className="m-0 font-mono text-xs">
          {node.mutations.length === 0 ? "none on this branch" : node.mutations.join(" ")}
        </dd>
      </dl>
      {rows.length > 0 && (
        <div className="border-line border-t px-3.5 py-3">
          <p className="mb-1.5 text-sm font-bold">
            Same {isTip ? "sample" : "clade"} in {rows.length} {rows.length === 1 ? "run" : "runs"}
          </p>
          <div className="plate rounded-md p-1.5">
            <DateIntervals rows={rows} onOpen={open} />
          </div>
          {missing > 0 && (
            <p className="text-ink-faint mt-1.5 text-xs">
              Not present in {missing} other time tree {missing === 1 ? "run" : "runs"}, whose trees lack this set of
              samples.
            </p>
          )}
        </div>
      )}
    </Panel>
  );
}

function useOtherTimetrees(runId: string): Array<{ id: string; label: string; clades: CladeIndex }> {
  const { data } = useRunList();

  const runs = useMemo(
    () => (data?.runs ?? []).filter((run) => run.id !== runId && run.command === "timetree" && run.status === "ok"),
    [data, runId],
  );

  const trees = useRunTrees(runs.map((run) => run.id));

  return runs.flatMap((run, index) => {
    const loaded = trees[index]?.data;

    return loaded === null || loaded === undefined ? [] : [{ id: run.id, label: run.title, clades: loaded.clades }];
  });
}
