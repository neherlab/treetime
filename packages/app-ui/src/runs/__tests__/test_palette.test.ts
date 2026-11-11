import { bottomAxisLabel, leftAxisLabel, niceAxis, niceStartAxis, yearTick } from "../palette";

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

  test("the domain of a nice axis covers values beyond the last round tick", () => {
    const axis = niceAxis([0, 1e-4, 3.2e-4]);

    expect(axis.domain[0]).toBeLessThanOrEqual(0);
    expect(axis.domain[1]).toBeGreaterThanOrEqual(3.2e-4);
  });

  test.each([
    [
      "a 49.25-year span takes the 1-2-5 step of 10 years and keeps no tick past the last value",
      [1964.4, 1990.2, 2013.65],
      { domain: [1960, 2013.65], ticks: [1960, 1970, 1980, 1990, 2000, 2010] },
    ],
    [
      "a last value on a round tick keeps that tick",
      [2000.1, 2014],
      { domain: [2000, 2014], ticks: [2000, 2002, 2004, 2006, 2008, 2010, 2012, 2014] },
    ],
  ])("a nice-start axis starts at a round tick and ends at the last value: %s", (_case, values, expected) => {
    expect(niceStartAxis(values)).toStrictEqual(expected);
  });
});

describe("chart axis labels", () => {
  test("a bottom axis label sits below the axis in the tick style", () => {
    expect(bottomAxisLabel("Date")).toStrictEqual({
      fontSize: 11,
      className: "fill-muted-foreground",
      value: "Date",
      position: "bottom",
      offset: 4,
    });
  });

  test("a left axis label reads upward inside the left edge in the tick style", () => {
    expect(leftAxisLabel("Branches")).toStrictEqual({
      fontSize: 11,
      className: "fill-muted-foreground",
      value: "Branches",
      angle: -90,
      position: "insideLeft",
    });
  });
});
