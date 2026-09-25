import { zAncestorShift } from "@neherlab/app-contracts";
import { useMemo } from "react";
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
} from "recharts";

import { formatSignedDays } from "../format";
import type { AncestorShift } from "../results/types";
import { PLATE, TICK_STYLE } from "./palette";

const HEIGHT = 300;

const MARGIN = { top: 8, right: 16, bottom: 24, left: 16 };

const DATA_EXTENT = ["dataMin", "dataMax"];

const CLADE_SIZE: [number, number] = [16, 160];

export function ShiftPlot({ shifts, firstLabel }: { shifts: readonly AncestorShift[]; firstLabel: string }) {
  const data = useMemo(
    () => shifts.map((shift) => ({ ...shift, x: shift.date_first.year, y: shift.shift_days, z: shift.tips })),
    [shifts],
  );

  return (
    <ResponsiveContainer width="100%" height={HEIGHT}>
      <ScatterChart margin={MARGIN}>
        <CartesianGrid stroke={PLATE.grid} />
        <XAxis type="number" dataKey="x" domain={DATA_EXTENT} tick={TICK_STYLE} tickFormatter={yearTick}>
          <Label value={`Date in ${firstLabel}`} position="bottom" offset={4} {...TICK_STYLE} />
        </XAxis>
        <YAxis type="number" dataKey="y" tick={TICK_STYLE} width={56}>
          <Label value="Shift in days" angle={-90} position="insideLeft" {...TICK_STYLE} />
        </YAxis>
        <ZAxis type="number" dataKey="z" range={CLADE_SIZE} />
        <ReferenceLine y={0} stroke={PLATE.muted} strokeDasharray="3 3" />
        <Tooltip content={<ShiftTooltip />} isAnimationActive={false} />
        <Scatter data={data} fill={PLATE.accent} fillOpacity={0.65} isAnimationActive={false} />
      </ScatterChart>
    </ResponsiveContainer>
  );
}

function ShiftTooltip({ active, payload }: { active?: boolean; payload?: ReadonlyArray<{ payload?: unknown }> }) {
  const parsed = zAncestorShift.safeParse(payload?.[0]?.payload);

  if (active !== true || !parsed.success) {
    return null;
  }

  const shift = parsed.data;

  return (
    <div className="rounded-md border border-[#bfcbc7] bg-white px-2.5 py-1.5 text-xs text-[#16302b] shadow-sm">
      <div className="font-bold">
        {shift.name}, {shift.tips} samples
      </div>
      <div>Date {shift.date_first.date}</div>
      <div>Shift {formatSignedDays(shift.shift_days)}</div>
    </div>
  );
}

function yearTick(value: number): string {
  return value.toFixed(1);
}
