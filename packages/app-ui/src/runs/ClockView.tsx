import type { RegressionStatistic, RootToTip, RunRecord, YearDate } from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { CellText, DataTable, dataColumns } from "@neherlab/design/ui/data-table";
import { DescriptionDetails, DescriptionList, DescriptionTerm } from "@neherlab/design/ui/description-list";
import { HoverCard, HoverCardContent, HoverCardTrigger } from "@neherlab/design/ui/hover-card";
import { useCallback, useMemo } from "react";

import { Panel } from "../components/Panel";
import { formatSignedDays, yearDateText } from "../format";
import { statisticText } from "../results/estimates";
import { OutputDownload } from "./OutputFiles";
import { clockSampleLists, useRootToTip } from "./rootToTip";
import { RootToTipPlot } from "./RootToTipPlot";
import { initialColorBy, TreeView, type TreeData } from "./TreeView";
import type { TreeLink } from "./TreeWorkspace";

interface SampleRow {
  name: string;
  date: YearDate;
  predictedDate: YearDate;
  residualDays: number | undefined;
  outlier: boolean;
}

const CLOCK_COLORINGS = ["num_date"];

const RESIDUAL_SORT = [{ id: "residual", desc: true }];

const sampleColumn = dataColumns<SampleRow>();

const SAMPLE_COLUMNS = sampleColumn.columns([
  sampleColumn.accessor((row) => row.name, {
    id: "name",
    header: "Sample",
    meta: { width: "2fr", minWidth: 240 },
    cell: ({ getValue }) => <CellText text={getValue()} />,
  }),
  sampleColumn.accessor((row) => row.date.year, {
    id: "date",
    header: "Date",
    meta: { width: 110, numeric: true },
    cell: ({ row }) => yearDateText(row.original.date),
  }),
  sampleColumn.accessor((row) => row.predictedDate.year, {
    id: "predicted",
    header: "Clock prediction",
    meta: { width: 130, numeric: true },
    cell: ({ row }) => yearDateText(row.original.predictedDate),
  }),
  sampleColumn.accessor((row) => (row.residualDays === undefined ? undefined : Math.abs(row.residualDays)), {
    id: "residual",
    header: "Residual",
    meta: { width: 100, numeric: true },
    sortUndefined: "last",
    cell: ({ row }) => (row.original.residualDays === undefined ? "" : formatSignedDays(row.original.residualDays)),
  }),
  sampleColumn.accessor((row) => (row.outlier ? "outlier" : "kept"), {
    id: "outlier",
    header: "Clock filter",
    meta: { width: 110 },
  }),
]);

export function ClockView({
  record,
  regression,
  statistics,
  tree,
}: {
  record: RunRecord;
  regression: RootToTip | undefined;
  statistics: readonly RegressionStatistic[];
  tree: TreeData | undefined;
}) {
  const { points, line } = useRootToTip(regression);

  const plot = useCallback(
    (link: TreeLink | undefined) => (
      <Panel
        figure
        title="Root-to-tip regression"
        hint="Red points are clock-filter outliers, grey points have an inferred date"
      >
        {points.length === 0 ? (
          <p className="text-muted-foreground px-2 py-3 text-sm">This run wrote no clock regression table.</p>
        ) : (
          <RootToTipPlot
            points={points}
            line={line}
            selected={link?.focus?.name}
            inView={link?.zoomed === true ? link.inView : undefined}
            onSelect={link?.select}
          />
        )}
      </Panel>
    ),
    [line, points],
  );

  const details = useCallback(
    (link: TreeLink | undefined) => (
      <>
        <RegressionStatistics record={record} statistics={statistics} />
        <ClockDetails regression={regression} onSelect={link?.select} />
      </>
    ),
    [record, regression, statistics],
  );

  return (
    <div className="grid gap-4">
      {tree === undefined ? (
        <>
          {plot(undefined)}
          {details(undefined)}
        </>
      ) : (
        <TreeView
          data={tree}
          colorBy={initialColorBy(tree.tree, CLOCK_COLORINGS)}
          primary="aside"
          aside={plot}
          below={details}
        />
      )}
    </div>
  );
}

function RegressionStatistics({
  record,
  statistics,
}: {
  record: RunRecord;
  statistics: readonly RegressionStatistic[];
}) {
  return (
    <Panel
      title="Regression statistics"
      hint="Root-to-tip regression of the clock model; hover a name for its definition"
      actions={
        <OutputDownload runId={record.id} kind="rtt">
          Root-to-tip table
        </OutputDownload>
      }
    >
      {statistics.length === 0 ? (
        <p className="text-muted-foreground px-3.5 py-3 text-sm">This run wrote no clock model.</p>
      ) : (
        <DescriptionList className="px-3.5 py-3">
          {statistics.map((statistic) => (
            <div key={statistic.statistic} className="contents">
              <DescriptionTerm>
                <HoverCard>
                  <HoverCardTrigger
                    render={
                      <button
                        type="button"
                        className="decoration-muted-foreground/60 cursor-help underline decoration-dotted underline-offset-4"
                      >
                        {statistic.label}
                      </button>
                    }
                  />
                  <HoverCardContent side="right" align="start" className="w-80 text-xs">
                    {statistic.definition}
                  </HoverCardContent>
                </HoverCard>
              </DescriptionTerm>
              <DescriptionDetails className="tabular-nums">{statisticText(statistic)}</DescriptionDetails>
            </div>
          ))}
        </DescriptionList>
      )}
    </Panel>
  );
}

function ClockDetails({
  regression,
  onSelect,
}: {
  regression: RootToTip | undefined;
  onSelect: ((name: string) => void) | undefined;
}) {
  const samples = useMemo(() => sampleRows(regression), [regression]);
  const lists = useMemo(() => clockSampleLists(regression), [regression]);

  return (
    <div className="grid gap-4">
      <div className="grid gap-4 @3xl:grid-cols-2">
        <SampleList
          title="Clock outliers"
          hint="Samples the clock filter left out of the clock model"
          names={lists.outliers}
          onSelect={onSelect}
        />
        <SampleList
          title="Undated samples"
          hint="Samples without a usable date; TreeTime inferred their dates"
          names={lists.undated}
          onSelect={onSelect}
        />
      </div>
      <Panel
        title="Residuals"
        hint="Residual = sampling date minus the date the clock model predicts from the root-to-tip distance"
      >
        <DataTable
          label="Residuals"
          columns={SAMPLE_COLUMNS}
          rows={samples}
          rowId={sampleKey}
          initialSorting={RESIDUAL_SORT}
          rowClassName={sampleTone}
        />
      </Panel>
    </div>
  );
}

function SampleList({
  title,
  hint,
  names,
  onSelect,
}: {
  title: string;
  hint: string;
  names: readonly string[];
  onSelect: ((name: string) => void) | undefined;
}) {
  return (
    <Panel title={`${title} (${names.length})`} hint={hint}>
      {names.length === 0 ? (
        <p className="text-muted-foreground px-3.5 py-3 text-sm">None</p>
      ) : (
        <ul className="max-h-60 overflow-y-auto px-2 py-2">
          {names.map((name) => (
            <li key={name}>
              <SampleButton name={name} onSelect={onSelect} />
            </li>
          ))}
        </ul>
      )}
    </Panel>
  );
}

function SampleButton({ name, onSelect }: { name: string; onSelect: ((name: string) => void) | undefined }) {
  const onClick = useCallback(() => onSelect?.(name), [name, onSelect]);

  if (onSelect === undefined) {
    return <span className="block px-1.5 py-0.5 font-mono text-xs">{name}</span>;
  }

  return (
    <Button
      type="button"
      variant="ghost"
      size="xs"
      className="h-auto w-full justify-start px-1.5 py-0.5 font-mono font-normal"
      onClick={onClick}
      title="Show in the tree and the plot"
    >
      {name}
    </Button>
  );
}

function sampleKey(row: SampleRow): string {
  return row.name;
}

function sampleTone(row: SampleRow): string | undefined {
  return row.outlier ? "bg-warning/10" : undefined;
}

function sampleRows(regression: RootToTip | undefined): SampleRow[] {
  return (regression?.points ?? []).flatMap((point) =>
    point.date === undefined
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
