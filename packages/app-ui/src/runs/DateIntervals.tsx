import { useCallback, useMemo } from "react";
import { Bar, ComposedChart, Scatter, XAxis, YAxis } from "recharts";

import type { DateInterval, YearDate } from "../results/types";
import { ChartTooltip, ChartTooltipFrame } from "../ui/chart";
import { CHART, niceAxis, THINNED_TICKS, TICK_STYLE, yearTick } from "./palette";

export interface DateRow {
  id: string;
  label: string;
  date: YearDate;
  interval: DateInterval | undefined;
  current: boolean;
}

const ROW_HEIGHT = 28;

const AXIS_HEIGHT = 28;

const MIN_SPAN_YEARS = 0.05;

const LABEL_WIDTH = 160;

const MARGIN = { top: 4, right: 16, bottom: 0, left: 0 };

const BAR_SIZE = 6;

const FULL_WIDTH = "100%";

export function DateIntervals({
  rows,
  onOpen,
}: {
  rows: readonly DateRow[];
  onOpen?: ((id: string) => void) | undefined;
}) {
  const axis = useMemo(() => dateAxis(rows), [rows]);
  const data = useMemo(() => rows.map(plotRow), [rows]);
  const labels = useMemo(() => new Map(rows.map((row) => [row.id, row.label])), [rows]);
  const rowLabel = useCallback((id: string) => labels.get(id) ?? id, [labels]);

  const openRow = useCallback(
    (point: { value?: unknown }) => {
      const row = rows.find((candidate) => candidate.id === point.value && !candidate.current);

      if (row !== undefined) {
        onOpen?.(row.id);
      }
    },
    [onOpen, rows],
  );

  return (
    <div className="[&_.recharts-cartesian-axis-tick_text]:fill-muted-foreground [&_.recharts-surface]:outline-hidden">
      <ComposedChart
        responsive
        layout="vertical"
        data={data}
        margin={MARGIN}
        width={FULL_WIDTH}
        height={rows.length * ROW_HEIGHT + AXIS_HEIGHT}
      >
        <XAxis
          type="number"
          domain={axis.domain}
          ticks={axis.ticks}
          {...THINNED_TICKS}
          tick={TICK_STYLE}
          tickFormatter={yearTick}
          stroke={CHART.muted}
        />
        <YAxis
          type="category"
          dataKey="id"
          width={LABEL_WIDTH}
          tick={TICK_STYLE}
          tickFormatter={rowLabel}
          onClick={openRow}
          className={onOpen === undefined ? "" : "cursor-pointer"}
        />
        <ChartTooltip content={<IntervalTooltip rows={rows} />} isAnimationActive={false} cursor={false} />
        <Bar dataKey="interval" barSize={BAR_SIZE} fill={CHART.accent} fillOpacity={0.35} isAnimationActive={false} />
        <Scatter dataKey="otherYear" fill={CHART.accent} isAnimationActive={false} />
        <Scatter dataKey="currentYear" fill={CHART.selection} isAnimationActive={false} />
      </ComposedChart>
    </div>
  );
}

function IntervalTooltip({
  rows,
  active,
  label,
}: {
  rows: readonly DateRow[];
  active?: boolean;
  label?: string | number;
}) {
  const row = rows.find((candidate) => candidate.id === label);

  if (active !== true || row === undefined) {
    return null;
  }

  return (
    <ChartTooltipFrame>
      <div className="font-medium">{row.label}</div>
      <div>{row.date.date}</div>
      {row.interval !== undefined && (
        <div className="text-muted-foreground">
          {row.interval.lower.date} to {row.interval.upper.date}
        </div>
      )}
    </ChartTooltipFrame>
  );
}

function plotRow(row: DateRow): PlotRow {
  return {
    id: row.id,
    otherYear: row.current ? null : row.date.year,
    currentYear: row.current ? row.date.year : null,
    interval: [row.interval?.lower.year ?? row.date.year, row.interval?.upper.year ?? row.date.year],
  };
}

function dateAxis(rows: readonly DateRow[]) {
  const ends = rows.flatMap((row) => [
    row.interval?.lower.year ?? row.date.year,
    row.interval?.upper.year ?? row.date.year,
  ]);

  const low = Math.min(...ends);
  const high = Math.max(...ends);
  const pad = Math.max((high - low) * 0.05, MIN_SPAN_YEARS / 2);

  return niceAxis([low - pad, high + pad]);
}

interface PlotRow {
  id: string;
  otherYear: number | null;
  currentYear: number | null;
  interval: [number, number];
}
