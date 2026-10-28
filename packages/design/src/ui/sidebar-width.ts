export const SIDEBAR_WIDTH_DEFAULT = 340;

export const SIDEBAR_WIDTH_MIN = 240;

export const SIDEBAR_WIDTH_MAX = 640;

export const SIDEBAR_WIDTH_STEP = 10;

export const MAIN_PANEL_MIN_WIDTH = 576;

const KEY_WIDTHS = new Map([
  ["ArrowLeft", (width: number) => width - SIDEBAR_WIDTH_STEP],
  ["ArrowRight", (width: number) => width + SIDEBAR_WIDTH_STEP],
  ["Home", () => SIDEBAR_WIDTH_MIN],
  ["End", () => SIDEBAR_WIDTH_MAX],
]);

export function clampSidebarWidth(width: number): number {
  return Math.round(Math.min(Math.max(width, SIDEBAR_WIDTH_MIN), SIDEBAR_WIDTH_MAX));
}

export function sidebarWidthOrDefault(width: number | undefined): number {
  return width === undefined ? SIDEBAR_WIDTH_DEFAULT : clampSidebarWidth(width);
}

export function fitSidebarWidth(width: number, viewportWidth: number): number {
  return viewportWidth > 0 ? clampSidebarWidth(Math.min(width, viewportWidth - MAIN_PANEL_MIN_WIDTH)) : width;
}

export function sidebarWidthForKey(width: number, key: string): number | undefined {
  const next = KEY_WIDTHS.get(key);

  return next === undefined ? undefined : clampSidebarWidth(next(width));
}
