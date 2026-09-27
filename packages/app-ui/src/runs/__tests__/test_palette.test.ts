import { describe, expect, test } from "vitest";

import { niceAxis, yearTick } from "../palette";

describe("chart axes", () => {
  test("year ticks drop trailing zeros", () => {
    expect([2016, 2016.5, 2013.25, 2017.8].map(yearTick)).toStrictEqual(["2016", "2016.5", "2013.25", "2017.8"]);
  });

  test("a nice axis covers the values with round ticks at both ends", () => {
    const axis = niceAxis([2017.77, 2024.6]);

    expect(axis.domain).toStrictEqual([axis.ticks.at(0), axis.ticks.at(-1)]);
    expect(axis.domain[0]).toBeLessThanOrEqual(2017.77);
    expect(axis.domain[1]).toBeGreaterThanOrEqual(2024.6);
    expect(axis.ticks.every(Number.isInteger)).toBe(true);
  });

  test("a nice axis extends only to the next round step beyond the values", () => {
    expect(niceAxis([2000.1, 2013.4])).toStrictEqual({
      domain: [2000, 2014],
      ticks: [2000, 2002, 2004, 2006, 2008, 2010, 2012, 2014],
    });
  });

  test("the domain of a nice axis covers values beyond the last round tick", () => {
    const axis = niceAxis([0, 1e-4, 3.2e-4]);

    expect(axis.domain[0]).toBeLessThanOrEqual(0);
    expect(axis.domain[1]).toBeGreaterThanOrEqual(3.2e-4);
  });
});
