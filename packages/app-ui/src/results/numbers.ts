import type { JsonFloat } from "@neherlab/app-contracts";

export function fromJsonFloat(value: JsonFloat): number {
  if (value === "inf") {
    return Number.POSITIVE_INFINITY;
  }

  if (value === "-inf") {
    return Number.NEGATIVE_INFINITY;
  }

  return value === "nan" ? Number.NaN : value;
}

export function nonFiniteLabel(value: number): string {
  if (Number.isNaN(value)) {
    return "not a number";
  }

  return value > 0 ? "+infinity" : "-infinity";
}
