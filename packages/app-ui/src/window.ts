import { MAIN_PANEL_MIN_WIDTH, SIDEBAR_WIDTH_MIN } from "./ui/sidebar-width";

export const WINDOW_BACKGROUND = { light: "#e9eeec", dark: "#0e1615" } as const;

export const WINDOW_MIN_SIZE = { width: SIDEBAR_WIDTH_MIN + MAIN_PANEL_MIN_WIDTH, height: 540 } as const;

export function windowBackground(dark: boolean): string {
  return dark ? WINDOW_BACKGROUND.dark : WINDOW_BACKGROUND.light;
}
