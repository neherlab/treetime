import type { ClockResults, RunRecord } from "@neherlab/app-contracts";

import { ClockView } from "./ClockView";
import type { TreeData } from "./TreeView";

export function ClockResultsView({
  record,
  data,
  tree,
}: {
  record: RunRecord;
  data: ClockResults;
  tree: TreeData | undefined;
}) {
  return <ClockView record={record} regression={data.root_to_tip} statistics={data.statistics} tree={tree} />;
}
