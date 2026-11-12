import type { RunRecord, TreeChange, TreeStatistics } from "@neherlab/app-contracts";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@neherlab/design/ui/table";

import { Panel } from "../components/Panel";
import { formatDuration } from "../format";
import { formatSignificant } from "../results/estimates";
import { MissingTree, TreeView, type TreeData } from "./TreeView";

const LENGTH_DIGITS = 4;

const STATISTICS: ReadonlyArray<{
  label: string;
  value: (statistics: TreeStatistics) => number | undefined;
  format: (value: number) => string;
}> = [
  { label: "Samples", value: (statistics) => statistics.samples, format: String },
  { label: "Internal nodes", value: (statistics) => statistics.internal_nodes, format: String },
  { label: "Polytomies", value: (statistics) => statistics.polytomies, format: String },
  {
    label: "Total branch length",
    value: (statistics) => statistics.total_branch_length,
    format: (value) => formatSignificant(value, LENGTH_DIGITS),
  },
];

export function TreeOnlyResults({
  record,
  data,
  tree,
}: {
  record: RunRecord;
  data: TreeChange;
  tree: TreeData | undefined;
}) {
  return (
    <div className="grid gap-4">
      {tree === undefined ? <MissingTree /> : <TreeView data={tree} colorBy={undefined} entropy={false} />}
      <Panel title="Input tree and result tree" hint={changeHint(record, data)}>
        <Table>
          <TableHeader>
            <TableRow>
              <TableHead>Statistic</TableHead>
              <TableHead className="text-right">Input tree</TableHead>
              <TableHead className="text-right">Result tree</TableHead>
              <TableHead className="text-right">Change</TableHead>
            </TableRow>
          </TableHeader>
          <TableBody>
            {STATISTICS.map((statistic) => {
              const before = data.before === undefined ? undefined : statistic.value(data.before);
              const after = data.after === undefined ? undefined : statistic.value(data.after);

              return (
                <TableRow key={statistic.label}>
                  <TableCell>{statistic.label}</TableCell>
                  <TableCell className="text-right tabular-nums">{valueText(before, statistic.format)}</TableCell>
                  <TableCell className="text-right tabular-nums">{valueText(after, statistic.format)}</TableCell>
                  <TableCell className="text-right tabular-nums">
                    {changeText(before, after, statistic.format)}
                  </TableCell>
                </TableRow>
              );
            })}
          </TableBody>
        </Table>
      </Panel>
    </div>
  );
}

function valueText(value: number | undefined, format: (value: number) => string): string {
  return value === undefined ? "-" : format(value);
}

function changeText(before: number | undefined, after: number | undefined, format: (value: number) => string): string {
  if (before === undefined || after === undefined) {
    return "-";
  }

  const change = after - before;

  return change === 0 ? "none" : `${change > 0 ? "+" : "-"}${format(Math.abs(change))}`;
}

function changeHint(record: RunRecord, data: TreeChange): string {
  const mutations = `${data.mutations} mutations on all branches of the result tree`;

  return record.duration_seconds === undefined
    ? mutations
    : `${mutations}; run time ${formatDuration(record.duration_seconds)}`;
}
