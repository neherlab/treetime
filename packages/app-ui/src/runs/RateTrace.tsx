import { useMemo } from "react";
import { Line, LineChart, ResponsiveContainer, Tooltip, XAxis, YAxis } from "recharts";

import type { IterationPoint } from "../results/progress";
import { CHART, TICK_STYLE } from "./palette";

const HEIGHT = 110;

const MARGIN = { top: 6, right: 12, bottom: 4, left: 8 };

const AUTO_EXTENT = ["auto", "auto"];

const DOT = { r: 2 };

export function RateTrace({ iterations }: { iterations: readonly IterationPoint[] }) {
  const data = useMemo(() => iterations.filter((point) => Number.isFinite(point.clockRate)), [iterations]);

  return (
    <ResponsiveContainer width="100%" height={HEIGHT}>
      <LineChart data={data} margin={MARGIN}>
        <XAxis dataKey="iteration" tick={TICK_STYLE} allowDecimals={false} />
        <YAxis domain={AUTO_EXTENT} tick={TICK_STYLE} tickFormatter={rateTick} width={56} />
        <Tooltip isAnimationActive={false} />
        <Line
          dataKey="clockRate"
          name="Clock rate"
          stroke={CHART.accent}
          strokeWidth={1.6}
          dot={DOT}
          isAnimationActive={false}
        />
      </LineChart>
    </ResponsiveContainer>
  );
}

function rateTick(value: number): string {
  return value.toExponential(2);
}
