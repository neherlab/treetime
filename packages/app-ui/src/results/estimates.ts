import type { RegressionStatistic } from "@neherlab/app-contracts";

import { yearDateText } from "../format";

const PLAIN_EXPONENTS = { lowest: -3, end: 4 } as const;

const DEFAULT_DIGITS = 3;

const SPREAD_DIGITS = 2;

const VALUE_DIGITS = { fewest: 2, most: 6 } as const;

const STATISTIC_DIGITS = 3;

export const RATE_UNIT = "per site per year";

export interface ValueWithSpread {
  value: string;
  spread: string | undefined;
}

export type Notation = "auto" | "plain" | "scientific";

export function formatSignificant(value: number, digits: number, notation: Notation = "auto"): string {
  if (value === 0 || !Number.isFinite(value)) {
    return String(value);
  }

  if (notation === "scientific" || (notation === "auto" && !isPlain(value))) {
    return value.toExponential(digits - 1).replace("e+", "e");
  }

  return value.toFixed(Math.max(0, digits - 1 - decimalExponent(value)));
}

export function valueWithSpread(
  value: number,
  spread: number | undefined,
  notation: Notation = "auto",
): ValueWithSpread {
  if (spread === undefined || !Number.isFinite(spread) || spread <= 0 || value === 0 || !Number.isFinite(value)) {
    return { value: formatSignificant(value, DEFAULT_DIGITS, notation), spread: undefined };
  }

  const digits = Math.min(
    VALUE_DIGITS.most,
    Math.max(VALUE_DIGITS.fewest, decimalExponent(value) - decimalExponent(spread) + 1),
  );

  const shared: Notation = notation === "auto" ? (isPlain(value) ? "plain" : "scientific") : notation;

  return {
    value: formatSignificant(value, digits, shared),
    spread: `± ${formatSignificant(spread, SPREAD_DIGITS, shared)}`,
  };
}

export function statisticText(statistic: RegressionStatistic): string {
  const value = statistic.value;

  if (value.kind === "date") {
    return yearDateText(value.date);
  }

  const number = formatSignificant(value.value, STATISTIC_DIGITS);

  return statistic.unit === undefined ? number : `${number} ${statistic.unit}`;
}

function isPlain(value: number): boolean {
  const exponent = decimalExponent(value);

  return exponent >= PLAIN_EXPONENTS.lowest && exponent < PLAIN_EXPONENTS.end;
}

function decimalExponent(value: number): number {
  return Math.floor(Math.log10(Math.abs(value)));
}
