import { describe, expect, test } from "vitest";

import { firstWindowBounds } from "../window-bounds";

const MIN_SIZE = { width: 816, height: 540 };

describe("window_bounds first start", () => {
  test("a laptop work area of 1280x672 gets 90% of it, centered", () => {
    expect(firstWindowBounds({ x: 0, y: 0, width: 1280, height: 672 }, MIN_SIZE)).toStrictEqual({
      x: 64,
      y: 33,
      width: 1152,
      height: 605,
    });
  });

  test("a large work area of 2560x1400 caps the window at 1600x1000", () => {
    expect(firstWindowBounds({ x: 0, y: 0, width: 2560, height: 1400 }, MIN_SIZE)).toStrictEqual({
      x: 480,
      y: 200,
      width: 1600,
      height: 1000,
    });
  });

  test("90% of a small work area below the minimum size is raised to the minimum size", () => {
    expect(firstWindowBounds({ x: 0, y: 0, width: 880, height: 580 }, MIN_SIZE)).toStrictEqual({
      x: 32,
      y: 20,
      width: 816,
      height: 540,
    });
  });

  test("a work area smaller than the minimum size limits the window to the work area", () => {
    expect(firstWindowBounds({ x: 0, y: 0, width: 800, height: 500 }, MIN_SIZE)).toStrictEqual({
      x: 0,
      y: 0,
      width: 800,
      height: 500,
    });
  });

  test("a work area offset from the origin keeps the window inside it", () => {
    expect(firstWindowBounds({ x: 1920, y: 25, width: 1440, height: 875 }, MIN_SIZE)).toStrictEqual({
      x: 1920 + 72,
      y: 25 + 43,
      width: 1296,
      height: 788,
    });
  });
});
