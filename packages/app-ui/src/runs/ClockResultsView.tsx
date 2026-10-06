import type { RunRecord, RunResults, ClockResults, YearDate } from "@neherlab/app-contracts";
import { useCallback, useMemo } from "react";

import { DataTable, dataColumns } from "../components/DataTable";
import { Panel, runTimeEntry, SummaryStrip, type SummaryEntry } from "../components/Panel";
import { formatRate, formatSignedDays, rSquaredText } from "../format";
import { OutputFiles } from "./OutputFiles";
import { useRootToTip } from "./rootToTip";
import { RootToTipPlot } from "./RootToTipPlot";
import { initialColorBy, MissingTree, TreeView, type TreeData } from "./TreeView";
import type { TreeLink } from "./TreeWorkspace";

interface SampleRow {
  name: string;
  date: YearDate;
  predictedDate: YearDate;
  residualDays: number;
  outlier: boolean;
}

const CLOCK_COLORINGS = ["num_date"];

const RESIDUAL_SORT = [{ id: "residual", desc: true }];

const sampleColumn = dataColumns<SampleRow>();

const SAMPLE_COLUMNS = sampleColumn.columns([
  sampleColumn.accessor((row) => row.name, { id: "name", header: "Sample", meta: { width: "1fr", minWidth: 160 } }),
  sampleColumn.accessor((row) => row.date.year, {
    id: "date",
    header: "Date",
    meta: { width: 110, numeric: true },
    cell: ({ row }) => row.original.date.date,
  }),
  sampleColumn.accessor((row) => row.predictedDate.year, {
    id: "predicted",
    header: "Clock prediction",
    meta: { width: 130, numeric: true },
    cell: ({ row }) => row.original.predictedDate.date,
  }),
  sampleColumn.accessor((row) => Math.abs(row.residualDays), {
    id: "residual",
    header: "Residual",
    meta: { width: 100, numeric: true },
    cell: ({ row }) => formatSignedDays(row.original.residualDays),
  }),
  sampleColumn.accessor((row) => (row.outlier ? "outlier" : "kept"), {
    id: "outlier",
    header: "Clock filter",
    meta: { width: 110 },
  }),
]);

export function ClockResultsView({
  record,
  results,
  data,
  tree,
}: {
  record: RunRecord;
  results: RunResults;
  data: ClockResults;
  tree: TreeData | undefined;
}) {
  const regression = data.root_to_tip;
  const samples = useMemo(() => sampleRows(data), [data]);

  const { points, line } = useRootToTip(regression);

  const summary = useMemo(() => clockSummary(record, data), [data, record]);

  const aside = useCallback(
    (link: TreeLink) => (
      <Panel figure title="Root-to-tip regression" hint="Dated samples; red points are clock-filter outliers">
        <RootToTipPlot
          points={points}
          line={line}
          selected={link.focus?.name}
          inView={undefined}
          onSelect={link.select}
        />
      </Panel>
    ),
    [line, points],
  );

  return (
    <div className="grid gap-4">
      <SummaryStrip entries={summary} />
      {tree === undefined ? (
        <MissingTree />
      ) : (
        <TreeView data={tree} colorBy={initialColorBy(tree.tree, CLOCK_COLORINGS)} aside={aside} />
      )}
      <Panel
        title="Samples"
        hint="Residual = sampling date minus the date the clock model predicts from the root-to-tip distance"
      >
        <DataTable
          label="Samples"
          columns={SAMPLE_COLUMNS}
          rows={samples}
          rowId={sampleKey}
          initialSorting={RESIDUAL_SORT}
          rowClassName={sampleTone}
        />
      </Panel>
      <OutputFiles record={record} citation={results.citation} />
    </div>
  );
}

function sampleKey(row: SampleRow): string {
  return row.name;
}

function sampleTone(row: SampleRow): string | undefined {
  return row.outlier ? "bg-warning/10" : undefined;
}

function sampleRows(data: ClockResults): SampleRow[] {
  return (data.root_to_tip?.points ?? []).flatMap((point) =>
    point.date === undefined || point.residual_days === undefined
      ? []
      : [
          {
            name: point.name,
            date: point.date,
            predictedDate: point.predicted_date,
            residualDays: point.residual_days,
            outlier: point.outlier,
          },
        ],
  );
}

function clockSummary(record: RunRecord, data: ClockResults): SummaryEntry[] {
  const estimates = data.estimates;
  const rate = estimates.clock_rate ?? undefined;
  const r = estimates.r ?? undefined;
  const outliers = estimates.outliers;

  return [
    {
      label: "Clock rate",
      value: rate === undefined ? "not written" : `${formatRate(rate)} /site/yr`,
      detail: estimates.clock_rate_fixed ? "Fixed" : "Root-to-tip regression",
    },
    {
      label: "Temporal signal",
      value: r === undefined ? "not computed" : `r = ${r.toFixed(3)}`,
      detail: rSquaredText(estimates.r_squared),
    },
    { label: "Dated samples", value: String(estimates.dated_samples) },
    {
      label: "Clock outliers",
      value: String(outliers),
      detail: outliers === 0 ? "None flagged by the clock filter" : "Flagged by the clock filter",
      tone: outliers === 0 ? undefined : "caution",
    },
    runTimeEntry(record),
  ];
}
