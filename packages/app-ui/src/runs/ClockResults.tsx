import type { RunRecord, RunResults } from "@neherlab/app-contracts";
import { useCallback, useMemo } from "react";

import { formatRate, formatSignedDays, rSquaredText } from "../format";
import type { ClockData, YearDate } from "../results/types";
import { OutputFiles } from "./OutputFiles";
import { Panel, runTimeEntry, SummaryStrip, type SummaryEntry } from "./Panel";
import { Plate } from "./Plate";
import { useRootToTip } from "./rootToTip";
import { RootToTipPlot } from "./RootToTipPlot";
import { SortableTable, type Column } from "./SortableTable";
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

const RESIDUAL_SORT = { key: "residual", descending: true };

const SAMPLE_COLUMNS: ReadonlyArray<Column<SampleRow>> = [
  { key: "name", label: "Sample", kind: "text", value: (row) => row.name },
  {
    key: "date",
    label: "Date",
    kind: "number",
    value: (row) => row.date.year,
    render: (row) => row.date.date,
  },
  {
    key: "predicted",
    label: "Clock prediction",
    kind: "number",
    value: (row) => row.predictedDate.year,
    render: (row) => row.predictedDate.date,
  },
  {
    key: "residual",
    label: "Residual",
    kind: "number",
    value: (row) => Math.abs(row.residualDays),
    render: (row) => formatSignedDays(row.residualDays),
  },
  { key: "outlier", label: "Clock filter", kind: "text", value: (row) => (row.outlier ? "outlier" : "kept") },
];

export function ClockResults({
  record,
  results,
  data,
  tree,
}: {
  record: RunRecord;
  results: RunResults;
  data: ClockData;
  tree: TreeData | undefined;
}) {
  const regression = data.root_to_tip;
  const samples = useMemo(() => sampleRows(data), [data]);

  const { points, line } = useRootToTip(regression);

  const summary = useMemo(() => clockSummary(record, data), [data, record]);

  const aside = useCallback(
    (link: TreeLink) => (
      <Plate title="Root-to-tip regression" caption="Dated samples; red points are clock-filter outliers">
        <RootToTipPlot
          points={points}
          line={line}
          selected={link.focus?.name}
          inView={undefined}
          onSelect={link.select}
        />
      </Plate>
    ),
    [line, points],
  );

  return (
    <div className="grid gap-3.5">
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
        <SortableTable
          label="Samples"
          columns={SAMPLE_COLUMNS}
          rows={samples}
          rowKey={sampleKey}
          initialSort={RESIDUAL_SORT}
          rowTone={sampleTone}
        />
      </Panel>
      <OutputFiles record={record} citation={results.citation} />
    </div>
  );
}

function sampleKey(row: SampleRow): string {
  return row.name;
}

function sampleTone(row: SampleRow): "caution" | undefined {
  return row.outlier ? "caution" : undefined;
}

function sampleRows(data: ClockData): SampleRow[] {
  return (data.root_to_tip?.points ?? []).flatMap((point) =>
    point.date === null || point.date === undefined || point.residual_days === null || point.residual_days === undefined
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

function clockSummary(record: RunRecord, data: ClockData): SummaryEntry[] {
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
