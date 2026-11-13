import type { RunRecord, TreeChange, TreeStatisticRow } from "@neherlab/app-contracts";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@neherlab/design/ui/table";

import { Panel } from "../components/Panel";
import { formatDuration } from "../format";
import { formatSignificant } from "../results/estimates";
import { MissingTree, TreeView, type TreeData } from "./TreeView";

const LENGTH_DIGITS = 4;

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
            {data.rows.map((row) => {
              const format = rowFormat(row);

              return (
                <TableRow key={row.label}>
                  <TableCell>{row.label}</TableCell>
                  <TableCell className="text-right tabular-nums">{valueText(row.before, format)}</TableCell>
                  <TableCell className="text-right tabular-nums">{valueText(row.after, format)}</TableCell>
                  <TableCell className="text-right tabular-nums">{changeText(row.before, row.after, format)}</TableCell>
                </TableRow>
              );
            })}
          </TableBody>
        </Table>
      </Panel>
    </div>
  );
}

function rowFormat(row: TreeStatisticRow): (value: number) => string {
  return row.count ? String : (value) => formatSignificant(value, LENGTH_DIGITS);
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
