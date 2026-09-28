import type { SettingSpec } from "./catalog";
import { cloneJson, getAt, sameJson, setAt, type JsonObject, type JsonValue } from "./json";

export function defaultConfig(specs: readonly SettingSpec[]): JsonObject {
  let config: JsonObject = {};

  for (const spec of specs) {
    config = setAt(config, spec.path, cloneJson(spec.default_value));
  }

  return config;
}

export function settingValue(config: JsonObject, spec: SettingSpec): JsonValue {
  return getAt(config, spec.path) ?? spec.default_value;
}

export function isChanged(config: JsonObject, spec: SettingSpec): boolean {
  return spec.role === "setting" && !sameJson(settingValue(config, spec), spec.default_value);
}

export function changedSpecs(specs: readonly SettingSpec[], config: JsonObject): SettingSpec[] {
  return specs.filter((spec) => isChanged(config, spec));
}

export function resetValue(spec: SettingSpec): JsonValue {
  return cloneJson(spec.default_value);
}

export function normalizeConfig(specs: readonly SettingSpec[], config: JsonObject): JsonObject {
  let normalized = defaultConfig(specs);

  for (const spec of specs) {
    const value = getAt(config, spec.path);

    if (value !== undefined && spec.role !== "output") {
      normalized = setAt(normalized, spec.path, cloneJson(value));
    }
  }

  return normalized;
}

export function carryOverConfig(
  specs: readonly SettingSpec[],
  previousSpecs: readonly SettingSpec[],
  previous: JsonObject,
): JsonObject {
  const carried = new Set(
    previousSpecs.flatMap((spec) => (spec.role === "input" || isChanged(previous, spec) ? [spec.key] : [])),
  );

  let config = defaultConfig(specs);

  for (const spec of specs) {
    const value = getAt(previous, spec.path);

    if (carried.has(spec.key) && value !== undefined) {
      config = setAt(config, spec.path, cloneJson(value));
    }
  }

  return config;
}
