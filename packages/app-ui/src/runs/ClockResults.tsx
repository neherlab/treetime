import type { RunRecordResult } from "@neherlab/app-contracts";
import { useCallback, useMemo } from "react";

import { daysBetween, formatDecimalDate, formatDuration, formatRate, formatSignedDays } from "../format";
import type { RunResults } from "../results/load";
import type { ClockRow } from "../results/readers";
import { initialColorBy } from "../results/tree";
import { OutputFiles } from "./OutputFiles";
import { Panel, SummaryStrip, type SummaryEntry } from "./Panel";
import { Plate } from "./Plate";
import { RootToTipPlot, type RttLine, type RttPoint } from "./RootToTipPlot";
import { SortableTable, type Column } from "./SortableTable";
import { MissingTree } from "./TimetreeResults";
import { TreeView } from "./TreeView";
import type { TreeLink } from "./TreeWorkspace";

interface SampleRow {
  name: string;
  div: number;
  date: number;
  predictedDate: number;
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
    value: (row) => row.date,
    render: (row) => formatDecimalDate(row.date),
  },
  {
    key: "predicted",
    label: "Clock prediction",
    kind: "number",
    value: (row) => row.predictedDate,
    render: (row) => formatDecimalDate(row.predictedDate),
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

export function ClockResults({ record, results }: { record: RunRecordResult; results: RunResults }) {
  const tips = useMemo(() => new Set(results.auspice?.tree.tips.map((tip) => tip.name) ?? []), [results.auspice]);
  const samples = useMemo(() => sampleRows(results.clockRows ?? [], tips), [results.clockRows, tips]);

  const points = useMemo<RttPoint[]>(
    () => samples.map((row) => ({ name: row.name, date: row.date, div: row.div, tip: true, excluded: row.outlier })),
    [samples],
  );

  const model = results.clockModel;

  const line = useMemo<RttLine | undefined>(
    () =>
      model === undefined
        ? undefined
        : {
            slope: model.rate,
            intercept: model.intercept,
            label: `TreeTime clock model: rate ${formatRate(model.rate)} /site/yr`,
          },
    [model],
  );

  const auspice = results.auspice;
  const summary = useMemo(() => clockSummary(record, results, samples), [record, results, samples]);

  const aside = useCallback(
    (link: TreeLink) => (
      <Plate title="Root-to-tip regression" caption="Dated samples; red rings mark clock-filter outliers">
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
      {auspice === undefined ? (
        <MissingTree />
      ) : (
        <TreeView auspice={auspice} colorBy={initialColorBy(auspice.tree, CLOCK_COLORINGS)} aside={aside} />
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
      <OutputFiles record={record} methods={undefined} />
    </div>
  );
}

function sampleKey(row: SampleRow): string {
  return row.name;
}

function sampleTone(row: SampleRow): "caution" | undefined {
  return row.outlier ? "caution" : undefined;
}

function sampleRows(rows: readonly ClockRow[], tips: ReadonlySet<string>): SampleRow[] {
  return rows.flatMap((row) =>
    row.date === undefined || !tips.has(row.name)
      ? []
      : [
          {
            name: row.name,
            div: row.div,
            date: row.date,
            predictedDate: row.predictedDate,
            residualDays: daysBetween(row.predictedDate, row.date),
            outlier: row.outlier,
          },
        ],
  );
}

function clockSummary(record: RunRecordResult, results: RunResults, samples: readonly SampleRow[]): SummaryEntry[] {
  const model = results.clockModel;
  const outliers = samples.filter((row) => row.outlier).length;

  return [
    {
      label: "Clock rate",
      value: model === undefined ? "not written" : `${formatRate(model.rate)} /site/yr`,
      detail: model?.fixed === true ? "Fixed" : "Root-to-tip regression",
    },
    {
      label: "Temporal signal",
      value: model?.r === undefined ? "not computed" : `r = ${model.r.toFixed(3)}`,
      detail: model?.r === undefined ? undefined : `R² = ${(model.r ** 2).toFixed(3)}`,
    },
    { label: "Dated samples", value: String(samples.length) },
    {
      label: "Clock outliers",
      value: String(outliers),
      detail: outliers === 0 ? "None flagged by the clock filter" : "Flagged by the clock filter",
      tone: outliers === 0 ? undefined : "caution",
    },
    {
      label: "Run time",
      value:
        record.duration_seconds === null || record.duration_seconds === undefined
          ? "-"
          : formatDuration(record.duration_seconds),
    },
  ];
}
