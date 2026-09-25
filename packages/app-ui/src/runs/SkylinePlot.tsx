import { useMemo } from "react";
import { Area, CartesianGrid, ComposedChart, Label, Line, ResponsiveContainer, Tooltip, XAxis, YAxis } from "recharts";

import type { SkylineSegment } from "../results/readers";
import { PLATE, TICK_STYLE } from "./palette";

const HEIGHT = 220;

const MARGIN = { top: 8, right: 16, bottom: 24, left: 16 };

const DATA_EXTENT = ["dataMin", "dataMax"];

const AUTO_EXTENT = ["auto", "auto"];

export function SkylinePlot({ segments }: { segments: readonly SkylineSegment[] }) {
  const data = useMemo(
    () =>
      segments.flatMap((segment) => [
        { t: segment.start, ne: segment.ne.value, band: [segment.ne.lower, segment.ne.upper] },
        { t: segment.end, ne: segment.ne.value, band: [segment.ne.lower, segment.ne.upper] },
      ]),
    [segments],
  );

  return (
    <ResponsiveContainer width="100%" height={HEIGHT}>
      <ComposedChart data={data} margin={MARGIN}>
        <CartesianGrid stroke={PLATE.grid} />
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
          fill={PLATE.accent}
          fillOpacity={0.15}
          isAnimationActive={false}
          name="Interval"
        />
        <Line
          dataKey="ne"
          type="linear"
          stroke={PLATE.accent}
          strokeWidth={1.8}
          dot={false}
          isAnimationActive={false}
          name="Ne"
        />
      </ComposedChart>
    </ResponsiveContainer>
  );
}

function yearTick(value: number): string {
  return value.toFixed(1);
}
