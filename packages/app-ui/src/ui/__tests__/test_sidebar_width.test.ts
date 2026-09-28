import { assert, double, integer, property, string } from "fast-check";
import { describe, expect, test } from "vitest";

import {
  clampSidebarWidth,
  parseSidebarWidth,
  SIDEBAR_WIDTH_DEFAULT,
  SIDEBAR_WIDTH_MAX,
  SIDEBAR_WIDTH_MIN,
  SIDEBAR_WIDTH_STEP,
  sidebarWidthForKey,
} from "../sidebar-width";

describe("sidebar width clamp", () => {
  test("a width inside the bounds is kept", () => {
    expect(clampSidebarWidth(300)).toBe(300);
  });

  test("a width below the minimum becomes the minimum", () => {
    expect(clampSidebarWidth(10)).toBe(SIDEBAR_WIDTH_MIN);
  });

  test("a width above the maximum becomes the maximum", () => {
    expect(clampSidebarWidth(5000)).toBe(SIDEBAR_WIDTH_MAX);
  });

  test("a fractional pointer position rounds to whole pixels", () => {
    expect(clampSidebarWidth(300.6)).toBe(301);
  });

  test("every finite width clamps to a whole pixel inside the bounds", () => {
    assert(
      property(double({ noNaN: true, noDefaultInfinity: true }), (width) => {
        const clamped = clampSidebarWidth(width);

        expect(Number.isInteger(clamped)).toBe(true);
        expect(clamped).toBeGreaterThanOrEqual(SIDEBAR_WIDTH_MIN);
        expect(clamped).toBeLessThanOrEqual(SIDEBAR_WIDTH_MAX);
      }),
    );
  });
});

describe("stored sidebar width", () => {
  test("a missing value gives the default width", () => {
    expect(parseSidebarWidth(undefined)).toBe(SIDEBAR_WIDTH_DEFAULT);
  });

  test("a stored number is restored", () => {
    expect(parseSidebarWidth("420")).toBe(420);
  });

  test("a non-numeric value gives the default width", () => {
    expect(parseSidebarWidth("wide")).toBe(SIDEBAR_WIDTH_DEFAULT);
  });

  test("a stored width outside the bounds is clamped", () => {
    expect(parseSidebarWidth("9000")).toBe(SIDEBAR_WIDTH_MAX);
  });

  test("any stored string restores a width inside the bounds", () => {
    assert(
      property(string(), (stored) => {
        const width = parseSidebarWidth(stored);

        expect(width).toBeGreaterThanOrEqual(SIDEBAR_WIDTH_MIN);
        expect(width).toBeLessThanOrEqual(SIDEBAR_WIDTH_MAX);
      }),
    );
  });
});

describe("sidebar width keys", () => {
  test("the default width is 4/3 of 16rem rounded to 10 pixels", () => {
    expect(SIDEBAR_WIDTH_DEFAULT).toBe(Math.round((256 * 4) / 3 / 10) * 10);
  });

  test("arrow right widens by one step", () => {
    expect(sidebarWidthForKey(300, "ArrowRight")).toBe(300 + SIDEBAR_WIDTH_STEP);
  });

  test("arrow left narrows by one step", () => {
    expect(sidebarWidthForKey(300, "ArrowLeft")).toBe(300 - SIDEBAR_WIDTH_STEP);
  });

  test("home and end jump to the bounds", () => {
    expect(sidebarWidthForKey(300, "Home")).toBe(SIDEBAR_WIDTH_MIN);
    expect(sidebarWidthForKey(300, "End")).toBe(SIDEBAR_WIDTH_MAX);
  });

  test("a step past a bound stops at the bound", () => {
    expect(sidebarWidthForKey(SIDEBAR_WIDTH_MAX, "ArrowRight")).toBe(SIDEBAR_WIDTH_MAX);
    expect(sidebarWidthForKey(SIDEBAR_WIDTH_MIN, "ArrowLeft")).toBe(SIDEBAR_WIDTH_MIN);
  });

  test("other keys leave the width unchanged", () => {
    expect(sidebarWidthForKey(300, "Enter")).toBeUndefined();
  });

  test("the arrow keys are inverse steps inside the bounds", () => {
    assert(
      property(integer({ min: SIDEBAR_WIDTH_MIN, max: SIDEBAR_WIDTH_MAX - SIDEBAR_WIDTH_STEP }), (width) => {
        const wider = sidebarWidthForKey(width, "ArrowRight") ?? Number.NaN;

        expect(sidebarWidthForKey(wider, "ArrowLeft")).toBe(width);
      }),
    );
  });
});
