import { useMemo } from "react";
import { Area, CartesianGrid, ComposedChart, Label, Line, ResponsiveContainer, Tooltip, XAxis, YAxis } from "recharts";

import type { SkylineSegment } from "../results/types";
import { CHART, PLOT_MARGIN, TICK_STYLE, yearTick } from "./palette";

const HEIGHT = 220;

const DATA_EXTENT = ["dataMin", "dataMax"];

const AUTO_EXTENT = ["auto", "auto"];

export function SkylinePlot({ segments }: { segments: readonly SkylineSegment[] }) {
  const data = useMemo(
    () =>
      segments.flatMap((segment) => [
        { t: segment.start, ne: segment.ne.value, band: band(segment) },
        { t: segment.end, ne: segment.ne.value, band: band(segment) },
      ]),
    [segments],
  );

  return (
    <ResponsiveContainer width="100%" height={HEIGHT}>
      <ComposedChart data={data} margin={PLOT_MARGIN}>
        <CartesianGrid stroke={CHART.grid} />
        <XAxis type="number" dataKey="t" domain={DATA_EXTENT} tick={TICK_STYLE} tickFormatter={yearTick}>
          <Label value="Date" position="bottom" offset={4} {...TICK_STYLE} />
        </XAxis>
        <YAxis type="number" scale="log" domain={AUTO_EXTENT} allowDataOverflow tick={TICK_STYLE} width={56}>
          <Label value="Effective population size" angle={-90} position="insideLeft" {...TICK_STYLE} />
        </YAxis>
        <Tooltip isAnimationActive={false} />
        <Area
          dataKey="band"
          type="linear"
          stroke="none"
          fill={CHART.accent}
          fillOpacity={0.15}
          isAnimationActive={false}
          name="Interval"
        />
        <Line
          dataKey="ne"
          type="linear"
          stroke={CHART.accent}
          strokeWidth={1.8}
          dot={false}
          isAnimationActive={false}
          name="Ne"
        />
      </ComposedChart>
    </ResponsiveContainer>
  );
}

function band(segment: SkylineSegment): [number, number] {
  return [segment.ne.lower ?? segment.ne.value, segment.ne.upper ?? segment.ne.value];
}
