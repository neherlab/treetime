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

export function getAt(object: SparseConfig, path: readonly string[]): JsonValue | undefined {
  let current: JsonValue | undefined = object;

  for (const key of path) {
    if (!isJsonObject(current)) {
      return undefined;
    }

    current = current[key];
  }

  return current;
}

export function setAt(object: SparseConfig, path: readonly string[], value: JsonValue | undefined): SparseConfig {
  const [head, ...rest] = path;

  if (head === undefined) {
    return object;
  }

  if (rest.length === 0) {
    if (value === undefined) {
      const { [head]: _removed, ...kept } = object;

      return kept;
    }

    return { ...object, [head]: value };
  }

  const child = object[head];

  return { ...object, [head]: setAt(isJsonObject(child) ? child : {}, rest, value) };
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
