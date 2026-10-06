import type { SparseConfig, JsonValue, SettingSpec } from "@neherlab/app-contracts";

import { cloneJson, getAt, sameJson, setAt } from "./json";

export function defaultConfig(specs: readonly SettingSpec[]): SparseConfig {
  let config: SparseConfig = {};

  for (const spec of specs) {
    if (spec.default_value !== undefined) {
      config = setAt(config, spec.path, cloneJson(spec.default_value));
    }
  }

  return config;
}

export function settingValue(config: SparseConfig, spec: SettingSpec): JsonValue | undefined {
  return getAt(config, spec.path) ?? spec.default_value;
}

export function isChanged(config: SparseConfig, spec: SettingSpec): boolean {
  return spec.role === "setting" && !sameJson(settingValue(config, spec), spec.default_value);
}

export function changedSpecs(specs: readonly SettingSpec[], config: SparseConfig): SettingSpec[] {
  return specs.filter((spec) => isChanged(config, spec));
}

export function resetValue(spec: SettingSpec): JsonValue | undefined {
  return spec.default_value === undefined ? undefined : cloneJson(spec.default_value);
}

export function normalizeConfig(specs: readonly SettingSpec[], config: SparseConfig): SparseConfig {
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
  previous: SparseConfig,
): SparseConfig {
  const carried = new Set(
    previousSpecs.flatMap((spec) => (spec.role === "input" || isChanged(previous, spec) ? [spec.key] : [])),
  );

  let config = defaultConfig(specs);

  for (const spec of specs) {
    const value = getAt(previous, spec.path);

    if (carried.has(spec.key) && value !== undefined && acceptsValue(spec, value)) {
      config = setAt(config, spec.path, cloneJson(value));
    }
  }

  return config;
}

function acceptsValue(spec: SettingSpec, value: JsonValue): boolean {
  if (spec.options.length === 0) {
    return true;
  }

  return (Array.isArray(value) ? value : [value]).every((item) =>
    spec.options.some((option) => sameJson(option.value, item)),
  );
}
