import type { JsonValue, SparseConfig } from "@neherlab/app-contracts";

export function isJsonObject(value: JsonValue | undefined): value is SparseConfig {
  return typeof value === "object" && !Array.isArray(value);
}

export function sameJson(left: JsonValue | undefined, right: JsonValue | undefined): boolean {
  if (left === undefined || right === undefined) {
    return left === right;
  }

  return canonicalJson(left) === canonicalJson(right);
}

function canonicalJson(value: JsonValue): string {
  return JSON.stringify(sortKeys(value));
}

export function withSetting(config: SparseConfig, key: string, value: JsonValue | undefined): SparseConfig {
  if (value === undefined) {
    const { [key]: _removed, ...kept } = config;

    return kept;
  }

  return { ...config, [key]: value };
}

export function cloneJson(value: JsonValue): JsonValue {
  return structuredClone(value);
}

function sortKeys(value: JsonValue): JsonValue {
  if (Array.isArray(value)) {
    return value.map(sortKeys);
  }

  if (isJsonObject(value)) {
    return Object.fromEntries(
      Object.entries(value)
        .toSorted(([left], [right]) => left.localeCompare(right))
        .map(([key, child]) => [key, sortKeys(child)]),
    );
  }

  return value;
}

export function isString(value: JsonValue | undefined): value is string {
  return typeof value === "string";
}

export function isNumber(value: JsonValue | undefined): value is number {
  return typeof value === "number";
}
