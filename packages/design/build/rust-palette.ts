import { rgbBytes } from "../src/tokens/color";
import { resolveColor, type ColorName } from "../src/tokens/tokens";

const CHART_ROLES = {
  PAPER: "chart-paper",
  INK: "chart-ink",
  MUTED: "chart-muted",
  GRID: "chart-grid",
  AXIS: "chart-axis",
  SERIES: "chart-series",
  FIT: "chart-fit",
  FAULT: "chart-fault",
} as const satisfies Record<string, ColorName>;

export function rustPalette(): string {
  const constants = Object.entries(CHART_ROLES).map(([constant, role]) => {
    const { r, g, b } = rgbBytes(resolveColor(role).light);

    return `pub const ${constant}: RGBColor = RGBColor(${hexByte(r)}, ${hexByte(g)}, ${hexByte(b)});`;
  });

  return ["use plotters::style::RGBColor;", "", constants.join("\n\n"), ""].join("\n");
}

function hexByte(value: number): string {
  return `0x${value.toString(16).padStart(2, "0")}`;
}
