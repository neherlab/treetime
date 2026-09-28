import { describe, expect, test } from "vitest";

import { FOLLOW_THRESHOLD_PX, isScrolledToEnd } from "../logFollow";

describe("log follow position", () => {
  test("a view scrolled to the last pixel is at the end", () => {
    expect(isScrolledToEnd({ scrollTop: 600, clientHeight: 400, scrollHeight: 1000 })).toBe(true);
  });

  test("content shorter than the view is at the end", () => {
    expect(isScrolledToEnd({ scrollTop: 0, clientHeight: 400, scrollHeight: 250 })).toBe(true);
  });

  test("a view within the threshold of the end is at the end", () => {
    expect(isScrolledToEnd({ scrollTop: 600 - FOLLOW_THRESHOLD_PX, clientHeight: 400, scrollHeight: 1000 })).toBe(true);
  });

  test("a view one pixel beyond the threshold is not at the end", () => {
    expect(isScrolledToEnd({ scrollTop: 600 - FOLLOW_THRESHOLD_PX - 1, clientHeight: 400, scrollHeight: 1000 })).toBe(
      false,
    );
  });

  test("a view at the top of long content is not at the end", () => {
    expect(isScrolledToEnd({ scrollTop: 0, clientHeight: 400, scrollHeight: 5000 })).toBe(false);
  });
});
