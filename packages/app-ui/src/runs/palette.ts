import { scaleLinear } from "d3-scale";

export const CHART = {
  ink: "var(--foreground)",
  muted: "var(--muted-foreground)",
  grid: "var(--chart-grid)",
  accent: "var(--chart-1)",
  fault: "var(--destructive)",
  tip: "var(--chart-3)",
  selection: "var(--chart-2)",
} as const;

export const TICK_STYLE = { fontSize: 11 } as const;

export const THINNED_TICKS = { interval: "equidistantPreserveStart", minTickGap: 56 } as const;

export const PLOT_MARGIN = { top: 8, right: 16, bottom: 24, left: 16 } as const;

const TICK_COUNT = 6;

export function yearTick(value: number): string {
  return String(Number(value.toFixed(2)));
}

export function niceAxis(values: readonly number[], tickCount: number = TICK_COUNT): AxisFrame {
  const min = Math.min(...values);
  const max = Math.max(...values);
  const scale = scaleLinear().domain([min, max]).nice(tickCount);
  const [low = min, high = max] = scale.domain();

  return { domain: [low, high], ticks: scale.ticks(tickCount) };
}

export interface AxisFrame {
  domain: [number, number];
  ticks: number[];
}

export function clamp(value: number, low: number, high: number): number {
  return Math.min(high, Math.max(low, value));
}
