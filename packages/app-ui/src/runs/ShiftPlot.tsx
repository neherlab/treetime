import { useElementSize } from "@mantine/hooks";
import { zAncestorShift } from "@neherlab/app-contracts";
import { useMemo } from "react";
import { CartesianGrid, Label, ReferenceLine, Scatter, ScatterChart, XAxis, YAxis, ZAxis } from "recharts";

import { formatSignedDays } from "../format";
import type { AncestorShift } from "../results/types";
import { ChartContainer, ChartTooltip } from "../ui/chart";
import { CHART, CHART_CONFIG, niceAxis, PLOT_MARGIN, TICK_STYLE, tickCountFor, yearTick } from "./palette";
import { TooltipCard } from "./TooltipCard";

const X_PX_PER_TICK = 100;

const CLADE_SIZE: [number, number] = [16, 160];

export function ShiftPlot({ shifts, firstLabel }: { shifts: readonly AncestorShift[]; firstLabel: string }) {
  const data = useMemo(
    () => shifts.map((shift) => ({ ...shift, x: shift.date_first.year, y: shift.shift_days, z: shift.tips })),
    [shifts],
  );

  const { ref: chart, width } = useElementSize<HTMLDivElement>();
  const xTickCount = tickCountFor(width, X_PX_PER_TICK);

  const xAxis = useMemo(
    () =>
      niceAxis(
        data.map((point) => point.x),
        xTickCount,
      ),
    [data, xTickCount],
  );

  return (
    <div ref={chart}>
      <ChartContainer config={CHART_CONFIG} className="aspect-auto h-[300px] w-full">
        <ScatterChart margin={PLOT_MARGIN}>
          <CartesianGrid stroke={CHART.grid} />
          <XAxis
            type="number"
            dataKey="x"
            domain={xAxis.domain}
            ticks={xAxis.ticks}
            tick={TICK_STYLE}
            tickFormatter={yearTick}
          >
            <Label
              value={`Date in ${firstLabel}`}
              position="bottom"
              offset={4}
              {...TICK_STYLE}
              className="fill-muted-foreground"
            />
          </XAxis>
          <YAxis type="number" dataKey="y" tick={TICK_STYLE} width={56}>
            <Label
              value="Shift in days"
              angle={-90}
              position="insideLeft"
              {...TICK_STYLE}
              className="fill-muted-foreground"
            />
          </YAxis>
          <ZAxis type="number" dataKey="z" range={CLADE_SIZE} />
          <ReferenceLine y={0} stroke={CHART.muted} strokeDasharray="3 3" />
          <ChartTooltip content={<ShiftTooltip />} isAnimationActive={false} />
          <Scatter data={data} fill={CHART.accent} fillOpacity={0.65} isAnimationActive={false} />
        </ScatterChart>
      </ChartContainer>
    </div>
  );
}

function ShiftTooltip({ active, payload }: { active?: boolean; payload?: ReadonlyArray<{ payload?: unknown }> }) {
  const parsed = zAncestorShift.safeParse(payload?.[0]?.payload);

  if (active !== true || !parsed.success) {
    return null;
  }

  const shift = parsed.data;

  return (
    <TooltipCard>
      <div className="font-semibold">
        {shift.name}, {shift.tips} samples
      </div>
      <div>Date {shift.date_first.date}</div>
      <div>Shift {formatSignedDays(shift.shift_days)}</div>
    </TooltipCard>
  );
}
