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

export function yearTick(value: number): string {
  return value.toFixed(1);
}
