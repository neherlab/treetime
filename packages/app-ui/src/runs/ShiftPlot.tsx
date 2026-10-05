import { zAncestorShift, type AncestorShift } from "@neherlab/app-contracts";
import { useMemo } from "react";
import { CartesianGrid, Label, ReferenceLine, Scatter, ScatterChart, XAxis, YAxis, ZAxis } from "recharts";

import { formatSignedDays } from "../format";
import { ChartContainer, ChartTooltip, ChartTooltipFrame } from "../ui/chart";
import { CHART, niceAxis, PLOT_MARGIN, THINNED_TICKS, TICK_STYLE, yearTick } from "./palette";

const CLADE_SIZE: [number, number] = [16, 160];

export function ShiftPlot({ shifts, firstLabel }: { shifts: readonly AncestorShift[]; firstLabel: string }) {
  const data = useMemo(
    () => shifts.map((shift) => ({ ...shift, x: shift.date_first.year, y: shift.shift_days, z: shift.tips })),
    [shifts],
  );

  const xAxis = useMemo(() => niceAxis(data.map((point) => point.x)), [data]);

  return (
    <ChartContainer className="aspect-auto h-[300px] w-full">
      <ScatterChart margin={PLOT_MARGIN}>
        <CartesianGrid stroke={CHART.grid} />
        <XAxis
          type="number"
          dataKey="x"
          domain={xAxis.domain}
          ticks={xAxis.ticks}
          {...THINNED_TICKS}
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
  );
}

function ShiftTooltip({ active, payload }: { active?: boolean; payload?: ReadonlyArray<{ payload?: unknown }> }) {
  const parsed = zAncestorShift.safeParse(payload?.[0]?.payload);

  if (active !== true || !parsed.success) {
    return null;
  }

  const shift = parsed.data;

  return (
    <ChartTooltipFrame>
      <div className="font-bold">
        {shift.name}, {shift.tips} samples
      </div>
      <div>Date {shift.date_first.date}</div>
      <div>Shift {formatSignedDays(shift.shift_days)}</div>
    </ChartTooltipFrame>
  );
}
