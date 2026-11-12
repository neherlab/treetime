import type { RegressionStatistic } from "@neherlab/app-contracts";

import { formatSignificant, statisticText, valueWithSpread } from "../estimates";

describe("numbers with significant digits", () => {
  test.each([
    [7.6123e-4, 3, "7.61e-4"],
    [0.91234, 3, "0.912"],
    [0.001234, 3, "0.00123"],
    [12.345, 3, "12.3"],
    [1234.5, 3, "1235"],
    [12_345, 3, "1.23e4"],
    [-2.5e-5, 2, "-2.5e-5"],
    [0, 3, "0"],
  ])("%d with %i digits reads %s", (value, digits, expected) => {
    expect(formatSignificant(value, digits)).toBe(expected);
  });
});

describe("a value with its spread", () => {
  test.each([
    { value: 7.6123e-4, spread: 4.1e-5, expected: { value: "7.6e-4", spread: "± 4.1e-5" } },
    { value: 7.6123e-4, spread: 4.1e-6, expected: { value: "7.61e-4", spread: "± 4.1e-6" } },
    { value: 7.6123e-4, spread: 2e-3, expected: { value: "7.6e-4", spread: "± 0.0020" } },
    { value: 1.234_567_89e-3, spread: 1e-12, expected: { value: "0.00123457", spread: "± 1.0e-12" } },
    { value: 7.6123e-4, spread: undefined, expected: { value: "7.61e-4", spread: undefined } },
    { value: 7.6123e-4, spread: 0, expected: { value: "7.61e-4", spread: undefined } },
  ])("$value with spread $spread reads $expected.value $expected.spread", ({ value, spread, expected }) => {
    expect(valueWithSpread(value, spread)).toStrictEqual(expected);
  });
});

describe("regression statistics", () => {
  const statistic = (value: RegressionStatistic["value"]): RegressionStatistic => ({
    statistic: "clock-rate",
    label: "Clock rate",
    definition: "Slope of the regression line.",
    value,
  });

  const withUnit = (value: RegressionStatistic["value"], unit: string): RegressionStatistic => ({
    ...statistic(value),
    unit,
  });

  test("a number shows three significant digits and its unit", () => {
    expect(statisticText(withUnit({ kind: "number", value: 7.6123e-4 }, "per site per year"))).toBe(
      "7.61e-4 per site per year",
    );
  });

  test("a number without a unit shows only the number", () => {
    expect(statisticText(statistic({ kind: "number", value: 0.83456 }))).toBe("0.835");
  });

  test("a date shows its calendar day", () => {
    expect(statisticText(statistic({ kind: "date", date: { year: 2013.84, date: "2013-11-03" } }))).toBe("2013-11-03");
  });
});
