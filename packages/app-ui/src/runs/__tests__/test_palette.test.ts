import { describe, expect, test } from "vitest";

import { niceAxis, tickCountFor, yearTick } from "../palette";

describe("chart axes", () => {
  test.each([
    [2016, "2016"],
    [2016.5, "2016.5"],
    [2013.25, "2013.25"],
    [2017.8, "2017.8"],
  ])("the year tick %d reads %s, without trailing zeros", (year, text) => {
    expect(yearTick(year)).toBe(text);
  });

  test("a nice axis covers the values with round ticks at both ends", () => {
    const axis = niceAxis([2017.77, 2024.6]);

    expect(axis.domain).toStrictEqual([axis.ticks.at(0), axis.ticks.at(-1)]);
    expect(axis.domain[0]).toBeLessThanOrEqual(2017.77);
    expect(axis.domain[1]).toBeGreaterThanOrEqual(2024.6);
    expect(axis.ticks.every(Number.isInteger)).toBe(true);
  });

  test("a 13.3-year span over 6 ticks takes the 1-2-5 step of 2 years and extends only to the next multiple of it", () => {
    expect(niceAxis([2000.1, 2013.4])).toStrictEqual({
      domain: [2000, 2014],
      ticks: [2000, 2002, 2004, 2006, 2008, 2010, 2012, 2014],
    });
  });

  test.each([
    [1100, 11],
    [0, 2],
    [99, 2],
    [5000, 12],
  ])("an axis %d px wide at 100 px per tick gets %d ticks, between 2 and 12", (width, count) => {
    expect(tickCountFor(width, 100)).toBe(count);
  });

  test("the domain of a nice axis covers values beyond the last round tick", () => {
    const axis = niceAxis([0, 1e-4, 3.2e-4]);

    expect(axis.domain[0]).toBeLessThanOrEqual(0);
    expect(axis.domain[1]).toBeGreaterThanOrEqual(3.2e-4);
  });
});
