import { scaleLinear } from "d3-scale";
import type { LabelProps } from "recharts";

export const CHART = {
  ink: "var(--chart-ink)",
  muted: "var(--chart-muted)",
  grid: "var(--chart-grid)",
  series: "var(--chart-series)",
  fit: "var(--chart-fit)",
  fault: "var(--chart-fault)",
  selection: "var(--chart-selection)",
} as const;

export const TICK_STYLE = { fontSize: 11 } as const;

export const THINNED_TICKS = { interval: "equidistantPreserveStart", minTickGap: 56 } as const;

export const PLOT_MARGIN = { top: 8, right: 16, bottom: 24, left: 16 } as const;

const AXIS_LABEL = { ...TICK_STYLE, className: "fill-muted-foreground" } as const;

export function bottomAxisLabel(value: string): LabelProps {
  return { ...AXIS_LABEL, value, position: "bottom", offset: 4 };
}

export function leftAxisLabel(value: string): LabelProps {
  return { ...AXIS_LABEL, value, angle: -90, position: "insideLeft" };
}

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
