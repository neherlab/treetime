import { hexColor } from "../src/tokens/color";
import { DATA_COLORS, resolveColor, type ColorName } from "../src/tokens/tokens";

export function matplotlibStyle(): string {
  const cycle = DATA_COLORS.map((name) => `'${color(name)}'`).join(", ");

  return [
    "# neherlab figure style for matplotlib: the colors and the typeface of the neherlab design language.",
    "# Generated from packages/design/src/tokens/tokens.ts by `just gen design`; edit the token source instead.",
    "# Use: matplotlib.style.use('path/to/neherlab.mplstyle'). Install IBM Plex Sans for the typeface;",
    "# without it, matplotlib falls back to DejaVu Sans.",
    "",
    `axes.prop_cycle: cycler('color', [${cycle}])`,
    "",
    "font.family: sans-serif",
    "font.sans-serif: IBM Plex Sans, DejaVu Sans, sans-serif",
    "font.monospace: IBM Plex Mono, DejaVu Sans Mono, monospace",
    "font.size: 10",
    "mathtext.fontset: dejavusans",
    "",
    `figure.facecolor: ${color("chart-paper")}`,
    `figure.edgecolor: ${color("chart-paper")}`,
    `savefig.facecolor: ${color("chart-paper")}`,
    `savefig.edgecolor: ${color("chart-paper")}`,
    "savefig.bbox: tight",
    "figure.dpi: 100",
    "savefig.dpi: 300",
    "",
    `axes.facecolor: ${color("chart-paper")}`,
    `axes.edgecolor: ${color("chart-axis")}`,
    "axes.linewidth: 0.8",
    `axes.labelcolor: ${color("chart-muted")}`,
    "axes.labelsize: 10",
    `axes.titlecolor: ${color("chart-ink")}`,
    "axes.titlesize: 12",
    "axes.titleweight: bold",
    "axes.titlelocation: left",
    "axes.spines.top: False",
    "axes.spines.right: False",
    "axes.grid: True",
    "axes.axisbelow: True",
    `grid.color: ${color("chart-grid")}`,
    "grid.linewidth: 0.8",
    "",
    `xtick.color: ${color("chart-axis")}`,
    `xtick.labelcolor: ${color("chart-muted")}`,
    "xtick.labelsize: 9",
    `ytick.color: ${color("chart-axis")}`,
    `ytick.labelcolor: ${color("chart-muted")}`,
    "ytick.labelsize: 9",
    "",
    `text.color: ${color("chart-ink")}`,
    "legend.frameon: False",
    "legend.fontsize: 9",
    `legend.labelcolor: ${color("chart-ink")}`,
    "",
    "lines.linewidth: 1.5",
    "lines.markersize: 5",
    "scatter.edgecolors: face",
    "image.cmap: cividis",
    "",
  ].join("\n");
}

function color(name: ColorName): string {
  return hexColor(resolveColor(name).light).slice(1);
}
