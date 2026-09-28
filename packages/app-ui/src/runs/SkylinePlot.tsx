import { useElementSize } from "@mantine/hooks";
import { useMemo } from "react";
import { Area, CartesianGrid, ComposedChart, Label, Line, XAxis, YAxis } from "recharts";

import type { SkylineSegment } from "../results/types";
import { ChartContainer, ChartTooltip, ChartTooltipContent } from "../ui/chart";
import { CHART, CHART_CONFIG, niceAxis, PLOT_MARGIN, TICK_STYLE, tickCountFor, yearTick } from "./palette";

const X_PX_PER_TICK = 100;

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

  const { ref: chart, width } = useElementSize<HTMLDivElement>();
  const xTickCount = tickCountFor(width, X_PX_PER_TICK);

  const xAxis = useMemo(
    () =>
      niceAxis(
        data.map((point) => point.t),
        xTickCount,
      ),
    [data, xTickCount],
  );

  return (
    <div ref={chart}>
      <ChartContainer config={CHART_CONFIG} className="aspect-auto h-[220px] w-full">
        <ComposedChart data={data} margin={PLOT_MARGIN}>
          <CartesianGrid stroke={CHART.grid} />
          <XAxis
            type="number"
            dataKey="t"
            domain={xAxis.domain}
            ticks={xAxis.ticks}
            tick={TICK_STYLE}
            tickFormatter={yearTick}
          >
            <Label value="Date" position="bottom" offset={4} {...TICK_STYLE} className="fill-muted-foreground" />
          </XAxis>
          <YAxis type="number" scale="log" domain={AUTO_EXTENT} allowDataOverflow tick={TICK_STYLE} width={56}>
            <Label
              value="Effective population size"
              angle={-90}
              position="insideLeft"
              {...TICK_STYLE}
              className="fill-muted-foreground"
            />
          </YAxis>
          <ChartTooltip content={<ChartTooltipContent />} isAnimationActive={false} />
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
      </ChartContainer>
    </div>
  );
}

function band(segment: SkylineSegment): [number, number] {
  return [segment.ne.lower ?? segment.ne.value, segment.ne.upper ?? segment.ne.value];
}
