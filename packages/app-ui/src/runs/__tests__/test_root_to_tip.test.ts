import type { RootToTipPoint } from "@neherlab/app-contracts";

import { clockSampleLists } from "../rootToTip";

const point = (name: string, extra: Partial<RootToTipPoint>): RootToTipPoint => ({
  name,
  div: 0.01,
  predicted_date: { year: 2015 },
  outlier: false,
  date: { year: 2015 },
  date_source: "input",
  ...extra,
});

describe("clock sample lists", () => {
  test("outliers and undated samples come from the regression points", () => {
    const missing: RootToTipPoint = {
      name: "missing",
      div: 0.01,
      predicted_date: { year: 2015 },
      outlier: false,
      date_source: "missing",
    };

    const regression = {
      points: [
        point("kept", {}),
        point("outlier", { outlier: true }),
        point("inferred", { date_source: "inferred" }),
        missing,
      ],
    };

    expect(clockSampleLists(regression)).toStrictEqual({
      outliers: ["outlier"],
      undated: ["inferred", "missing"],
    });
  });

  test("a run without a regression table has empty lists", () => {
    expect(clockSampleLists(undefined)).toStrictEqual({ outliers: [], undated: [] });
  });
});
