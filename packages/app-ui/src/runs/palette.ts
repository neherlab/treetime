import { getNiceTickValues } from "recharts";

export const CHART = {
  ink: "var(--color-ink)",
  muted: "var(--color-ink-muted)",
  faint: "var(--color-ink-faint)",
  grid: "var(--color-chart-grid)",
  accent: "var(--color-accent)",
  fault: "var(--color-signal-danger)",
  tip: "var(--color-chart-tip)",
  selection: "var(--color-chart-selection)",
} as const;

export const TICK_STYLE = { fontSize: 11, fill: CHART.muted } as const;

export const PLOT_MARGIN = { top: 8, right: 16, bottom: 24, left: 16 } as const;

const TICK_COUNT = 6;

export function yearTick(value: number): string {
  return String(Number(value.toFixed(2)));
}

export function niceAxis(values: readonly number[], tickCount: number = TICK_COUNT): AxisFrame {
  const low = Math.min(...values);
  const high = Math.max(...values);
  const ticks = getNiceTickValues([low, high], tickCount, true);

  return { domain: [Math.min(ticks.at(0) ?? low, low), Math.max(ticks.at(-1) ?? high, high)], ticks };
}

export interface AxisFrame {
  domain: [number, number];
  ticks: number[];
}
