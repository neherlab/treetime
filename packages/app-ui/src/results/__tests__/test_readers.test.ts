import { describe, expect, test } from "vitest";

import { readClockModel, readTracelog } from "../readers";

describe("output readers", () => {
  test("a tracelog with non-finite likelihoods keeps them", () => {
    const text = [
      "n_diff,n_resolved,max_time_change,rms_time_change,log_lh_seq,log_lh_pos,log_lh_coal,log_lh_total",
      "0,0,0.1,0.01,-10,-2,inf,inf",
      "0,0,NaN,0.01,-10,-2,,-12",
    ].join("\n");

    expect(readTracelog(text).map((row) => [row.logLhCoal, row.logLhTotal, row.maxTimeChange])).toStrictEqual([
      [Number.POSITIVE_INFINITY, Number.POSITIVE_INFINITY, 0.1],
      [undefined, -12, Number.NaN],
    ]);
  });

  test("a fixed clock has no regression statistics", () => {
    expect(readClockModel('{"clock_rate": 0.001, "intercept": -2.0, "stats": "fixed"}')).toStrictEqual({
      rate: 0.001,
      intercept: -2,
      fixed: true,
      r: undefined,
    });
  });
});
