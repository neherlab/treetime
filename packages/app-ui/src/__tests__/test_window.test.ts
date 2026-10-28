import { MAIN_PANEL_MIN_WIDTH, SIDEBAR_WIDTH_MIN } from "@neherlab/design/ui/sidebar-width";

import { WINDOW_BACKGROUND, WINDOW_MIN_SIZE, windowBackground } from "../window";

const HEX_COLOR = /^#[0-9a-f]{6}$/u;

describe("window background", () => {
  test("the window colors are opaque hex colors, the format Electron accepts", () => {
    expect([WINDOW_BACKGROUND.light, WINDOW_BACKGROUND.dark]).toStrictEqual([
      expect.stringMatching(HEX_COLOR),
      expect.stringMatching(HEX_COLOR),
    ]);
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
