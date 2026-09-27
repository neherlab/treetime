import { useCallback, useId, useMemo, useState } from "react";
import {
  CartesianGrid,
  Label,
  ReferenceLine,
  ResponsiveContainer,
  Scatter,
  ScatterChart,
  Tooltip,
  XAxis,
  YAxis,
  ZAxis,
  type ScatterPointItem,
} from "recharts";
import * as z from "zod";

import { Switch } from "../ui";
import { ChartTooltip } from "./ChartTooltip";
import { CHART, PLOT_MARGIN, TICK_STYLE, yearTick } from "./palette";
import { type PlacedPoint, type PlotFrame, placePoints, plotFrame } from "./rttFrame";

export interface RttPoint {
  name: string;
  date: number;
  dateText: string;
  div: number;
  excluded: boolean;
  inferred: boolean;
}

export interface RttLine {
  slope: number;
  intercept: number;
  label: string;
}

const HEIGHT = 300;

const FADED_OPACITY = 0.18;

const OFF_AXES_OPACITY = 0.45;

const TIP_SIZE: [number, number] = [36, 36];

const RING_SIZE: [number, number] = [220, 220];

const SERIES_LOOK: readonly SeriesLook[] = [
  { key: "inferredFaded", fill: CHART.faint, opacity: FADED_OPACITY },
  { key: "tipsFaded", fill: CHART.tip, opacity: FADED_OPACITY },
  { key: "excludedFaded", fill: CHART.fault, opacity: FADED_OPACITY },
  { key: "inferred", fill: CHART.faint, opacity: 1 },
  { key: "tips", fill: CHART.tip, opacity: 1 },
  { key: "excluded", fill: CHART.fault, opacity: 1 },
];

const zPointPayload = z.object({
  name: z.string(),
  date: z.number(),
  dateText: z.string(),
  div: z.number(),
  excluded: z.boolean(),
  inferred: z.boolean(),
  offAxes: z.boolean(),
});

export function RootToTipPlot({
  points,
  line,
  selected,
  inView,
  onSelect,
}: {
  points: readonly RttPoint[];
  line: RttLine | undefined;
  selected: string | undefined;
  inView: ReadonlySet<string> | undefined;
  onSelect: ((name: string) => void) | undefined;
}) {
  const [fitToModel, setFitToModel] = useState(true);
  const fitId = useId();
  const hasOutliers = useMemo(() => points.some((point) => point.excluded), [points]);
  const frame = useMemo(() => plotFrame(points, fitToModel && hasOutliers), [fitToModel, hasOutliers, points]);
  const placed = useMemo(() => placePoints(points, frame), [frame, points]);
  const series = useMemo(() => pointSeries(placed, selected, inView), [inView, placed, selected]);
  const segment = useMemo(() => lineSegment(frame, line), [frame, line]);

  return (
    <div>
      <div className="flex flex-wrap items-center justify-between gap-x-4 gap-y-1 px-2 pb-1 text-xs">
        {line === undefined ? (
          <span />
        ) : (
          <p className="text-ink-muted m-0">
            <span className="bg-accent mr-1.5 inline-block h-0.5 w-4 align-middle" />
            {line.label}
          </p>
        )}
        {hasOutliers && (
          <label htmlFor={fitId} className="text-ink-muted flex cursor-pointer items-center gap-2">
            <Switch id={fitId} checked={fitToModel} onCheckedChange={setFitToModel} />
            Fit axes to the samples in the clock model
          </label>
        )}
      </div>
      <ResponsiveContainer width="100%" height={HEIGHT}>
        <ScatterChart margin={PLOT_MARGIN}>
          <CartesianGrid stroke={CHART.grid} />
          <XAxis
            type="number"
            dataKey="x"
            domain={frame.x.domain}
            ticks={frame.x.ticks}
            tick={TICK_STYLE}
            tickFormatter={yearTick}
            stroke={CHART.faint}
          >
            <Label value="Date" position="bottom" offset={4} {...TICK_STYLE} />
          </XAxis>
          <YAxis
            type="number"
            dataKey="y"
            domain={frame.y.domain}
            ticks={frame.y.ticks}
            tick={TICK_STYLE}
            tickFormatter={divergenceTick}
            stroke={CHART.faint}
            width={56}
          >
            <Label value="Divergence from the root" angle={-90} position="insideLeft" {...TICK_STYLE} />
          </YAxis>
          <ZAxis zAxisId="tip" range={TIP_SIZE} />
          <ZAxis zAxisId="ring" range={RING_SIZE} />
          <Tooltip content={<PointTooltip />} isAnimationActive={false} />
          {segment !== undefined && (
            <ReferenceLine segment={segment} stroke={CHART.accent} strokeWidth={1.5} ifOverflow="hidden" />
          )}
          {SERIES_LOOK.map((look) => (
            <SeriesScatter key={look.key} look={look} points={series[look.key]} onSelect={onSelect} />
          ))}
          <OffAxesScatter points={series.offAxes} onSelect={onSelect} />
          <Scatter
            data={series.selected}
            zAxisId="ring"
            fill="none"
            stroke={CHART.selection}
            strokeWidth={2}
            isAnimationActive={false}
          />
        </ScatterChart>
      </ResponsiveContainer>
    </div>
  );
}

function SeriesScatter({
  look,
  points,
  onSelect,
}: {
  look: SeriesLook;
  points: readonly PlacedPoint[];
  onSelect: ((name: string) => void) | undefined;
}) {
  const select = useSelect(points, onSelect);

  return (
    <Scatter
      data={points}
      zAxisId="tip"
      fill={look.fill}
      fillOpacity={look.opacity}
      isAnimationActive={false}
      onClick={select}
    />
  );
}

function OffAxesScatter({
  points,
  onSelect,
}: {
  points: readonly PlacedPoint[];
  onSelect: ((name: string) => void) | undefined;
}) {
  const select = useSelect(points, onSelect);

  return (
    <Scatter
      data={points}
      zAxisId="tip"
      fill={CHART.fault}
      fillOpacity={0}
      stroke={CHART.fault}
      strokeOpacity={OFF_AXES_OPACITY}
      strokeWidth={1.5}
      isAnimationActive={false}
      onClick={select}
    />
  );
}

function useSelect(points: readonly PlacedPoint[], onSelect: ((name: string) => void) | undefined) {
  return useCallback(
    (_item: ScatterPointItem, index: number) => {
      const point = points[index];

      if (point !== undefined) {
        onSelect?.(point.name);
      }
    },
    [onSelect, points],
  );
}

function pointSeries(
  points: readonly PlacedPoint[],
  selected: string | undefined,
  inView: ReadonlySet<string> | undefined,
) {
  const visible = (point: PlacedPoint) => inView === undefined || inView.has(point.name);
  const onAxes = points.filter((point) => !point.offAxes);

  const of = (kind: PointRole, shown: boolean) =>
    onAxes.filter((point) => pointRole(point) === kind && visible(point) === shown);

  return {
    inferred: of("inferred", true),
    inferredFaded: of("inferred", false),
    tips: of("tip", true),
    tipsFaded: of("tip", false),
    excluded: of("excluded", true),
    excludedFaded: of("excluded", false),
    offAxes: points.filter((point) => point.offAxes),
    selected: points.filter((point) => point.name === selected),
  };
}

function pointRole(point: RttPoint): PointRole {
  if (point.excluded) {
    return "excluded";
  }

  return point.inferred ? "inferred" : "tip";
}

type PointRole = "inferred" | "tip" | "excluded";

type SeriesKey = "inferred" | "inferredFaded" | "tips" | "tipsFaded" | "excluded" | "excludedFaded";

interface SeriesLook {
  key: SeriesKey;
  fill: string;
  opacity: number;
}

function PointTooltip({ active, payload }: { active?: boolean; payload?: ReadonlyArray<{ payload?: unknown }> }) {
  const point = zPointPayload.safeParse(payload?.[0]?.payload);

  if (active !== true || !point.success) {
    return null;
  }

  return (
    <ChartTooltip>
      <div className="font-bold">{point.data.name}</div>
      <div>Date {point.data.dateText}</div>
      <div>Divergence {point.data.div.toExponential(3)}</div>
      {point.data.inferred && <div>Date inferred by the time tree; the sample has no input date</div>}
      {point.data.excluded && <div className="text-signal-danger">Excluded from the clock model</div>}
      {point.data.offAxes && <div className="text-ink-muted">Outside the axes; drawn at their edge</div>}
    </ChartTooltip>
  );
}

function lineSegment(frame: PlotFrame, line: RttLine | undefined) {
  if (line === undefined) {
    return undefined;
  }

  const [from, to] = frame.x.domain;

  return [
    { x: from, y: line.slope * from + line.intercept },
    { x: to, y: line.slope * to + line.intercept },
  ] as const;
}

function divergenceTick(value: number): string {
  return value === 0 ? "0" : value.toExponential(1);
}
