import type { RunRecord, TimetreeEstimates, TimetreeResults } from "@neherlab/app-contracts";
import { useCallback } from "react";

import { Panel } from "../components/Panel";
import { useRootToTip } from "./rootToTip";
import { RootToTipPlot } from "./RootToTipPlot";
import { SelectionPanel } from "./SelectionPanel";
import { TimetreeEstimatesPanel } from "./TimetreeEstimates";
import { initialColorBy, MissingTree, TreeView, type TreeData } from "./TreeView";
import type { TreeLink } from "./TreeWorkspace";

const TIMETREE_COLORINGS = ["num_date"];

export function TimetreeResultsView({
  record,
  data,
  tree,
}: {
  record: RunRecord;
  data: TimetreeResults;
  tree: TreeData | undefined;
}) {
  const estimates = data.estimates;

  if (tree === undefined || estimates === undefined) {
    return <MissingTree />;
  }

  return <TimetreeView record={record} data={data} estimates={estimates} tree={tree} />;
}

function TimetreeView({
  record,
  data,
  estimates,
  tree,
}: {
  record: RunRecord;
  data: TimetreeResults;
  estimates: TimetreeEstimates;
  tree: TreeData;
}) {
  const aside = useCallback(
    (link: TreeLink) => <TimetreeAside record={record} tree={tree} data={data} link={link} />,
    [data, record, tree],
  );

  return (
    <div className="grid gap-4">
      <TimetreeEstimatesPanel record={record} estimates={estimates} />
      <TreeView data={tree} colorBy={initialColorBy(tree.tree, TIMETREE_COLORINGS)} entropy={false} aside={aside} />
    </div>
  );
}

function TimetreeAside({
  record,
  tree,
  data,
  link,
}: {
  record: RunRecord;
  tree: TreeData;
  data: TimetreeResults;
  link: TreeLink;
}) {
  const { points, line } = useRootToTip(data.root_to_tip);
  const root = tree.tree.nodes[0];
  const selected = link.focus !== undefined && link.focus !== root;

  return (
    <>
      {selected && <SelectionPanel runId={record.id} title={record.title} tree={tree.tree} link={link} />}
      <Panel
        figure
        title="Root-to-tip regression"
        hint="Samples as the final clock model saw them; red points are clock-filter outliers, grey points have an inferred date"
      >
        {points.length === 0 ? (
          <p className="text-muted-foreground px-2 py-3 text-sm">This run wrote no clock regression table.</p>
        ) : (
          <RootToTipPlot
            points={points}
            line={line}
            selected={link.focus?.name}
            inView={link.zoomed ? link.inView : undefined}
            onSelect={link.select}
          />
        )}
      </Panel>
    </>
  );
}
