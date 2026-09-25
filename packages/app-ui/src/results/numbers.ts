export type JsonFloat = number | "inf" | "-inf" | "nan";

const NON_FINITE = new Map([
  ["inf", Number.POSITIVE_INFINITY],
  ["+inf", Number.POSITIVE_INFINITY],
  ["infinity", Number.POSITIVE_INFINITY],
  ["-inf", Number.NEGATIVE_INFINITY],
  ["-infinity", Number.NEGATIVE_INFINITY],
  ["nan", Number.NaN],
]);

export function fromJsonFloat(value: JsonFloat): number {
  if (value === "inf") {
    return Number.POSITIVE_INFINITY;
  }

  if (value === "-inf") {
    return Number.NEGATIVE_INFINITY;
  }

  return value === "nan" ? Number.NaN : value;
}

export function parseNumber(text: string): number {
  const trimmed = text.trim();
  const special = NON_FINITE.get(trimmed.toLowerCase());

  if (special !== undefined) {
    return special;
  }

  const value = Number(trimmed);

  if (trimmed === "" || Number.isNaN(value)) {
    throw new Error(`"${text}" is not a number`);
  }

  return value;
}

export function parseOptionalNumber(text: string | undefined): number | undefined {
  return text === undefined || text.trim() === "" ? undefined : parseNumber(text);
}

export function nonFiniteLabel(value: number): string {
  if (Number.isNaN(value)) {
    return "not a number";
  }

  return value > 0 ? "+infinity" : "-infinity";
}
