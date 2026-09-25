import { useCallback, useMemo } from "react";
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

import { ChartTooltip } from "./ChartTooltip";
import { PLATE, PLOT_MARGIN, TICK_STYLE, yearTick } from "./palette";

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

const DATA_EXTENT = ["dataMin", "dataMax"];

const TIP_SIZE: [number, number] = [36, 36];

const RING_SIZE: [number, number] = [220, 220];

const TIP_COLOR = "#27477f";

const SERIES_LOOK: readonly SeriesLook[] = [
  { key: "inferredFaded", fill: PLATE.faint, opacity: FADED_OPACITY },
  { key: "tipsFaded", fill: TIP_COLOR, opacity: FADED_OPACITY },
  { key: "excludedFaded", fill: PLATE.fault, opacity: FADED_OPACITY },
  { key: "inferred", fill: PLATE.faint, opacity: 1 },
  { key: "tips", fill: TIP_COLOR, opacity: 1 },
  { key: "excluded", fill: PLATE.fault, opacity: 1 },
];

const zPointPayload = z.object({
  name: z.string(),
  date: z.number(),
  dateText: z.string(),
  div: z.number(),
  excluded: z.boolean(),
  inferred: z.boolean(),
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
  const series = useMemo(() => pointSeries(points, selected, inView), [inView, points, selected]);
  const segment = useMemo(() => lineSegment(points, line), [line, points]);

  return (
    <div>
      {line !== undefined && (
        <p className="text-plate-muted m-0 px-2 pb-1 text-xs">
          <span className="bg-plate-accent mr-1.5 inline-block h-0.5 w-4 align-middle" />
          {line.label}
        </p>
      )}
      <ResponsiveContainer width="100%" height={HEIGHT}>
        <ScatterChart margin={PLOT_MARGIN}>
          <CartesianGrid stroke={PLATE.grid} />
          <XAxis type="number" dataKey="date" domain={DATA_EXTENT} tick={TICK_STYLE} tickFormatter={yearTick}>
            <Label value="Date" position="bottom" offset={4} {...TICK_STYLE} />
          </XAxis>
          <YAxis type="number" dataKey="div" tick={TICK_STYLE} tickFormatter={divergenceTick} width={56}>
            <Label value="Divergence from the root" angle={-90} position="insideLeft" {...TICK_STYLE} />
          </YAxis>
          <ZAxis zAxisId="tip" range={TIP_SIZE} />
          <ZAxis zAxisId="ring" range={RING_SIZE} />
          <Tooltip content={<PointTooltip />} isAnimationActive={false} />
          {segment !== undefined && (
            <ReferenceLine segment={segment} stroke={PLATE.accent} strokeWidth={1.5} ifOverflow="extendDomain" />
          )}
          {SERIES_LOOK.map((look) => (
            <SeriesScatter key={look.key} look={look} points={series[look.key]} onSelect={onSelect} />
          ))}
          <Scatter
            data={series.selected}
            zAxisId="ring"
            fill="none"
            stroke={PLATE.selection}
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
  points: readonly RttPoint[];
  onSelect: ((name: string) => void) | undefined;
}) {
  const select = useCallback(
    (_item: ScatterPointItem, index: number) => {
      const point = points[index];

      if (point !== undefined) {
        onSelect?.(point.name);
      }
    },
    [onSelect, points],
  );

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

function pointSeries(
  points: readonly RttPoint[],
  selected: string | undefined,
  inView: ReadonlySet<string> | undefined,
) {
  const visible = (point: RttPoint) => inView === undefined || inView.has(point.name);

  const of = (kind: PointRole, shown: boolean) =>
    points.filter((point) => pointRole(point) === kind && visible(point) === shown);

  return {
    inferred: of("inferred", true),
    inferredFaded: of("inferred", false),
    tips: of("tip", true),
    tipsFaded: of("tip", false),
    excluded: of("excluded", true),
    excludedFaded: of("excluded", false),
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
      {point.data.excluded && <div className="text-plate-fault">Excluded from the clock model</div>}
    </ChartTooltip>
  );
}

function lineSegment(points: readonly RttPoint[], line: RttLine | undefined) {
  if (line === undefined || points.length === 0) {
    return undefined;
  }

  const dates = points.map((point) => point.date);
  const from = Math.min(...dates);
  const to = Math.max(...dates);

  return [
    { x: from, y: line.slope * from + line.intercept },
    { x: to, y: line.slope * to + line.intercept },
  ] as const;
}

function divergenceTick(value: number): string {
  return value === 0 ? "0" : value.toExponential(1);
}
