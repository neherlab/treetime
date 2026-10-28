import { hexPair } from "@neherlab/design/color";
import { BASE_COLORS } from "@neherlab/design/tokens";
import { MAIN_PANEL_MIN_WIDTH, SIDEBAR_WIDTH_MIN } from "@neherlab/design/ui/sidebar-width";

export const WINDOW_BACKGROUND = hexPair(BASE_COLORS.background);

export const WINDOW_MIN_SIZE = { width: SIDEBAR_WIDTH_MIN + MAIN_PANEL_MIN_WIDTH, height: 540 } as const;

export function windowBackground(dark: boolean): string {
  return dark ? WINDOW_BACKGROUND.dark : WINDOW_BACKGROUND.light;
}
