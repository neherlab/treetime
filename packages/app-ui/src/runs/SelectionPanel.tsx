import { useNavigate } from "@tanstack/react-router";
import { useCallback, useMemo } from "react";

import { formatLevel } from "../format";
import { useCladeInRuns } from "../queries";
import type { ResultTree } from "../results/types";
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
  const root = tree.nodes[0];
  const node = link.focus ?? root;
  const navigate = useNavigate();
  const { data: found } = useCladeInRuns(runId, node?.name ?? "");
  const open = useCallback((id: string) => void navigate({ to: "/runs/$id/results", params: { id } }), [navigate]);

  const rows = useMemo(() => {
    const current: DateRow[] =
      node?.date === null || node?.date === undefined
        ? []
        : [{ id: runId, label: title, date: node.date, interval: node.date_interval ?? undefined, current: true }];

    return [
      ...current,
      ...(found?.matches ?? []).flatMap((match) =>
        match.date === null || match.date === undefined
          ? []
          : [
              {
                id: match.run,
                label: match.title,
                date: match.date,
                interval: match.date_interval ?? undefined,
                current: false,
              },
            ],
      ),
    ];
  }, [found, node, runId, title]);

  if (node === undefined) {
    return null;
  }

  const missing = (found?.searched_runs ?? 0) - (found?.matches.length ?? 0);
  const interval = node.date_interval ?? undefined;
  const isTip = node.children.length === 0;

  return (
    <Panel
      title={isTip ? node.name : `Ancestor of ${node.tips} samples`}
      hint={
        link.focus === undefined || link.focus === root
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
        <dd className="m-0 tabular-nums">
          {node.date === null || node.date === undefined ? "not dated" : node.date.date}
        </dd>
        {interval !== undefined && (
          <>
            <dt className="text-ink-faint">{formatLevel(interval.level)} interval</dt>
            <dd className="m-0 tabular-nums">
              {interval.lower.date} to {interval.upper.date} ({Math.round(interval.days)} days)
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
