import type { SkylineSegment } from "@neherlab/app-contracts";
import { ChartContainer, ChartTooltip, ChartTooltipContent } from "@neherlab/design/ui/chart";
import { useMemo } from "react";
import { Area, CartesianGrid, ComposedChart, Line, XAxis, YAxis } from "recharts";

import {
  bottomAxisLabel,
  CHART,
  leftAxisLabel,
  niceAxis,
  PLOT_MARGIN,
  THINNED_TICKS,
  TICK_STYLE,
  yearTick,
} from "./palette";

const AUTO_EXTENT = ["auto", "auto"];

export function SkylinePlot({ segments }: { segments: readonly SkylineSegment[] }) {
  const data = useMemo(
    () =>
      segments.flatMap((segment) => [
        { t: segment.start, tc: segment.tc.value, band: band(segment) },
        { t: segment.end, tc: segment.tc.value, band: band(segment) },
      ]),
    [segments],
  );

  const xAxis = useMemo(() => niceAxis(data.map((point) => point.t)), [data]);

  return (
    <ChartContainer className="aspect-auto h-[220px] w-full">
      <ComposedChart data={data} margin={PLOT_MARGIN}>
        <CartesianGrid stroke={CHART.grid} />
        <XAxis
          type="number"
          dataKey="t"
          domain={xAxis.domain}
          ticks={xAxis.ticks}
          {...THINNED_TICKS}
          tick={TICK_STYLE}
          tickFormatter={yearTick}
          label={bottomAxisLabel("Date")}
        />
        <YAxis
          type="number"
          scale="log"
          domain={AUTO_EXTENT}
          allowDataOverflow
          tick={TICK_STYLE}
          width={56}
          label={leftAxisLabel("Tc (years)")}
        />
        <ChartTooltip content={<ChartTooltipContent />} isAnimationActive={false} />
        <Area
          dataKey="band"
          type="linear"
          stroke="none"
          fill={CHART.series}
          fillOpacity={0.15}
          isAnimationActive={false}
          name="Interval"
        />
        <Line
          dataKey="tc"
          type="linear"
          stroke={CHART.series}
          strokeWidth={1.8}
          dot={false}
          isAnimationActive={false}
          name="Tc"
        />
      </ComposedChart>
    </ChartContainer>
  );
}

function band(segment: SkylineSegment): [number, number] {
  return [segment.tc.lower ?? segment.tc.value, segment.tc.upper ?? segment.tc.value];
}
