import { describe, expect, test } from "vitest";

import { leastSquares } from "../regression";

describe("least squares line", () => {
  test("points on a line give that line", () => {
    expect(leastSquares([0, 1, 2, 3].map((x) => ({ x, y: 2 * x - 1 })))).toStrictEqual({ slope: 2, intercept: -1 });
  });

  test("symmetric scatter around a line gives the line", () => {
    expect(
      leastSquares([
        { x: 0, y: 1 },
        { x: 0, y: -1 },
        { x: 2, y: 3 },
        { x: 2, y: 1 },
      ]),
    ).toStrictEqual({ slope: 1, intercept: 0 });
  });

  test("a single point or points at one date give no line", () => {
    expect([
      leastSquares([{ x: 1, y: 1 }]),
      leastSquares([
        { x: 1, y: 1 },
        { x: 1, y: 2 },
      ]),
    ]).toStrictEqual([undefined, undefined]);
  });
});
