import { describe, expect, test } from "vitest";

import type { RttPoint } from "../RootToTipPlot";
import { placePoints, plotFrame } from "../rttFrame";

const KEPT: readonly RttPoint[] = [
  { name: "a", date: 2018.2, dateText: "2018-03-14", div: 1e-4, excluded: false, inferred: false },
  { name: "b", date: 2020.5, dateText: "2020-07-01", div: 2e-4, excluded: false, inferred: false },
  { name: "c", date: 2023.9, dateText: "2023-11-25", div: 3.2e-4, excluded: false, inferred: true },
];

const OUTLIER: RttPoint = {
  name: "far",
  date: 2021.4,
  dateText: "2021-05-27",
  div: 5.6e-3,
  excluded: true,
  inferred: false,
};

const POINTS = [...KEPT, OUTLIER];

describe("root-to-tip plot frame", () => {
  test("fitting to the clock model leaves the outlier out of the axis ranges", () => {
    const frame = plotFrame(POINTS, true);

    expect(frame.y.domain[1]).toBeLessThan(OUTLIER.div);
    expect(frame.y.domain[1]).toBeGreaterThanOrEqual(3.2e-4);
    expect(frame.x.domain[0]).toBeLessThanOrEqual(2018.2);
    expect(frame.x.domain[1]).toBeGreaterThanOrEqual(2023.9);
  });

  test("without fitting, the axis ranges cover every point", () => {
    const frame = plotFrame(POINTS, false);

    expect(frame.y.domain[1]).toBeGreaterThanOrEqual(OUTLIER.div);
  });

  test("the divergence axis starts at zero or below", () => {
    expect(plotFrame(POINTS, true).y.domain[0]).toBeLessThanOrEqual(0);
  });

  test.each(["x", "y"] as const)("the %s ticks start and end at the axis range", (key) => {
    const axis = plotFrame(POINTS, true)[key];

    expect([axis.ticks.at(0), axis.ticks.at(-1)]).toStrictEqual(axis.domain);
  });

  test("the year ticks fall on whole years for a multi-year range", () => {
    expect(plotFrame(POINTS, true).x.ticks.every(Number.isInteger)).toBe(true);
  });

  test("when every point is an outlier, the axes cover all points", () => {
    const frame = plotFrame([OUTLIER], true);

    expect(frame.y.domain[1]).toBeGreaterThanOrEqual(OUTLIER.div);
  });

  test("a point outside the axes is drawn at their edge and keeps its values", () => {
    const frame = plotFrame(POINTS, true);
    const placed = placePoints(POINTS, frame);
    const far = placed.find((point) => point.name === "far");

    expect(far).toMatchObject({ offAxes: true, x: OUTLIER.date, y: frame.y.domain[1], div: OUTLIER.div });
    expect(placed.filter((point) => point.offAxes).map((point) => point.name)).toStrictEqual(["far"]);
  });

  test("points inside the axes keep their coordinates", () => {
    const placed = placePoints(KEPT, plotFrame(KEPT, true));

    expect(placed.map((point) => [point.x, point.y, point.offAxes])).toStrictEqual(
      KEPT.map((point) => [point.date, point.div, false]),
    );
  });
});
