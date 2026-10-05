import { describe, expect, test } from "vitest";

import { MAIN_PANEL_MIN_WIDTH, SIDEBAR_WIDTH_MIN } from "../ui/sidebar-width";
import { WINDOW_BACKGROUND, WINDOW_MIN_SIZE, windowBackground } from "../window";

import themeCss from "../theme.css?raw";

describe("window background", () => {
  test("the window colors equal the page backgrounds of the light and the dark theme", () => {
    expect(WINDOW_BACKGROUND).toStrictEqual({
      light: themeBackground(":root,\n.light-scope {"),
      dark: themeBackground(".dark {"),
    });
  });

  test.each([
    [true, WINDOW_BACKGROUND.dark],
    [false, WINDOW_BACKGROUND.light],
  ])("a dark native theme %s picks %s", (dark, expected) => {
    expect(windowBackground(dark)).toBe(expected);
  });
});

describe("window minimum size", () => {
  test("the minimum width holds the narrowest sidebar next to the narrowest main panel", () => {
    expect(WINDOW_MIN_SIZE.width).toBe(SIDEBAR_WIDTH_MIN + MAIN_PANEL_MIN_WIDTH);
  });

  test("the minimum size fits the smallest common laptop work area of 1097x569", () => {
    expect([WINDOW_MIN_SIZE.width <= 1097, WINDOW_MIN_SIZE.height <= 569]).toStrictEqual([true, true]);
  });
});

function themeBackground(selector: string): string | undefined {
  const block = themeCss.slice(themeCss.indexOf(selector));

  return /--background:\s*(#[0-9a-f]{6});/u.exec(block)?.[1];
}
