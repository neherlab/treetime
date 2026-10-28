import {
  formatHex,
  modeLrgb,
  modeOklab,
  modeOklch,
  modeRgb,
  toGamut,
  useMode,
  wcagContrast,
  type Oklch as CuloriOklch,
  type Rgb,
} from "culori/fn";

import type { Oklch } from "./tokens";

const MODES = [modeRgb, modeLrgb, modeOklab, modeOklch] as const;

for (const mode of MODES) {
  // oxlint-disable-next-line react-hooks/rules-of-hooks -- culori's useMode registers a color space; it is not a React hook
  useMode(mode);
}

const toSrgb = toGamut("rgb", "oklch");

export function cssColor(color: Oklch): string {
  return `oklch(${trim(color.l)} ${trim(color.c)} ${trim(color.h)})`;
}

export function hexColor(color: Oklch): string {
  return formatHex(srgb(color));
}

export function rgbBytes(color: Oklch): RgbBytes {
  const { r, g, b } = srgb(color);

  return { r: byte(r), g: byte(g), b: byte(b) };
}

export function contrastRatio(foreground: Oklch, background: Oklch): number {
  return wcagContrast(srgb(foreground), srgb(background));
}

export interface RgbBytes {
  r: number;
  g: number;
  b: number;
}

function srgb(color: Oklch): Rgb {
  const culoriColor: CuloriOklch = { mode: "oklch", l: color.l, c: color.c, h: color.h };

  return toSrgb(culoriColor);
}

function byte(channel: number): number {
  return Math.round(Math.min(1, Math.max(0, channel)) * 255);
}

function trim(value: number): string {
  return String(Number(value.toFixed(4)));
}
