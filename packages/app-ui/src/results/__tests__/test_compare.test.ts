import { describe, expect, test } from "vitest";

import { indexClades, matchAncestors } from "../clades";
import { ancestorShifts, compareEstimates, timetreeEstimates, type TimetreeEstimates } from "../estimates";
import { parseAuspiceJson, readAuspiceTree, type ResultTree } from "../tree";

const FIRST = tree({
  meta: {},
  tree: node("root", 2000, [
    node("AB", 2001, [tip("A", 2003), tip("B", 2004)]),
    node("CD", 2001.5, [tip("C", 2003), tip("D", 2004)]),
  ]),
});

const SECOND = tree({
  meta: {},
  tree: node("top", 2000.5, [
    node("DC", 2001.25, [tip("D", 2004), tip("C", 2003)]),
    node("ABE", 2000.75, [tip("B", 2004), tip("A", 2003), tip("E", 2002)]),
  ]),
});

describe("clade matching", () => {
  test("a clade with its tips in another order matches", () => {
    const matched = matchAncestors(indexClades(FIRST), indexClades(SECOND));

    expect(matched.map((pair) => [pair.first.name, pair.second.name])).toStrictEqual([["CD", "DC"]]);
  });

  test("a clade present in one run only has no match", () => {
    const names = new Set(matchAncestors(indexClades(FIRST), indexClades(SECOND)).map((pair) => pair.first.name));

    expect({ ab: names.has("AB"), root: names.has("root") }).toStrictEqual({ ab: false, root: false });
  });
});

describe("run comparison", () => {
  test("shifts are the calendar days between the dates of matched ancestors", () => {
    expect(ancestorShifts(FIRST, SECOND)).toStrictEqual([
      { name: "CD", tips: 2, dateFirst: 2001.5, shiftDays: -91.25 },
    ]);
  });

  test("the root shift, interval widths and rate change are reported in days and percent", () => {
    const first = estimates({ rootDate: 2015, rootInterval: [2014.75, 2015.25], rate: 0.001 });
    const second = estimates({ rootDate: 2015.5, rootInterval: [2015, 2016], rate: 0.0011 });
    const comparison = compareEstimates(first, second);

    expect({
      rootShiftDays: comparison.rootShiftDays,
      intervalWidthDays: comparison.intervalWidthDays,
      ratePercentChange: comparison.ratePercentChange?.toFixed(6),
    }).toStrictEqual({ rootShiftDays: 182.5, intervalWidthDays: [182.5, 365], ratePercentChange: "10.000000" });
  });

  test("a missing estimate leaves its difference empty", () => {
    const comparison = compareEstimates(
      estimates({ rootDate: undefined, rootInterval: undefined, rate: undefined }),
      estimates({}),
    );

    expect(comparison).toStrictEqual({
      rootShiftDays: undefined,
      intervalWidthDays: [undefined, 182.5],
      ratePercentChange: undefined,
    });
  });
});

describe("timetree estimates", () => {
  test("a root date within 5% of the interval width from a bound is near the edge", () => {
    expect([2015.04, 2015.05, 2015.06, 2015.97].map((date) => rootNearEdge(date, [2015, 2016]))).toStrictEqual([
      true,
      true,
      false,
      true,
    ]);
  });

  test("a root without an interval is never near the edge", () => {
    expect(rootNearEdge(2015.5, [2015.5, 2015.5])).toBe(false);
  });
});

function estimates(overrides: Partial<TimetreeEstimates>): TimetreeEstimates {
  return {
    rootDate: 2015,
    rootInterval: [2015, 2015.5],
    rootNearIntervalEdge: false,
    rate: 0.001,
    rateStd: undefined,
    rateFixed: false,
    r: undefined,
    samples: 4,
    excludedSamples: 0,
    logLikelihood: undefined,
    iterations: 0,
    ...overrides,
  };
}

function rootNearEdge(date: number, confidence: [number, number]): boolean {
  const rooted = tree({
    meta: {},
    tree: { name: "root", node_attrs: { num_date: { value: date, confidence } }, children: [tip("A", 2016)] },
  });

  return timetreeEstimates({ tree: rooted, clockModel: undefined, augurClock: undefined, trace: undefined })
    .rootNearIntervalEdge;
}

function node(name: string, date: number, children: unknown[]) {
  return { name, node_attrs: { div: 0, num_date: { value: date } }, children };
}

function tip(name: string, date: number) {
  return { name, node_attrs: { div: 0, num_date: { value: date } } };
}

function tree(document: { meta: Record<string, never>; tree: unknown }): ResultTree {
  return readAuspiceTree(parseAuspiceJson(JSON.stringify(document)));
}
