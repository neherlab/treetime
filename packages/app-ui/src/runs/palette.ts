export const PLATE = {
  ink: "#16302b",
  muted: "#4b5f5a",
  faint: "#8a9a95",
  grid: "#e3e9e7",
  accent: "#17695a",
  fault: "#b42318",
  selection: "#882255",
} as const;

export const TICK_STYLE = { fontSize: 11, fill: PLATE.muted } as const;

export const PLOT_MARGIN = { top: 8, right: 16, bottom: 24, left: 16 } as const;

export function yearTick(value: number): string {
  return value.toFixed(1);
}
