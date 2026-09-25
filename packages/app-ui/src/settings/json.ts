import * as z from "zod";

export type JsonValue = string | number | boolean | null | JsonValue[] | JsonObject;

export interface JsonObject {
  [key: string]: JsonValue;
}

export const zJsonValue: z.ZodType<JsonValue> = z.lazy(() =>
  z.union([z.string(), z.number(), z.boolean(), z.null(), z.array(zJsonValue), zJsonObject]),
);

export const zJsonObject: z.ZodType<JsonObject> = z.lazy(() => z.record(z.string(), zJsonValue));

export function isJsonObject(value: JsonValue | undefined): value is JsonObject {
  return typeof value === "object" && value !== null && !Array.isArray(value);
}

export function sameJson(left: JsonValue | undefined, right: JsonValue | undefined): boolean {
  return canonicalJson(left ?? null) === canonicalJson(right ?? null);
}

function canonicalJson(value: JsonValue): string {
  return JSON.stringify(sortKeys(value));
}

export function getAt(object: JsonObject, path: readonly string[]): JsonValue | undefined {
  let current: JsonValue | undefined = object;

  for (const key of path) {
    if (!isJsonObject(current)) {
      return undefined;
    }

    current = current[key];
  }

  return current;
}

export function setAt(object: JsonObject, path: readonly string[], value: JsonValue): JsonObject {
  const [head, ...rest] = path;

  if (head === undefined) {
    return object;
  }

  if (rest.length === 0) {
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
      Object.keys(value)
        .toSorted()
        .map((key) => [key, sortKeys(value[key] ?? null)]),
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
