import { scaleLinear } from "d3-scale";

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

const MIN_TICKS = 2;

const MAX_TICKS = 12;

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

export function tickCountFor(width: number, pxPerTick: number): number {
  return clamp(Math.floor(width / pxPerTick), MIN_TICKS, MAX_TICKS);
}
