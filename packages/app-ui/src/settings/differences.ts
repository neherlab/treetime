import type { RunInput } from "@neherlab/app-contracts";

import type { SettingSpec } from "./catalog";
import { settingValue } from "./config";
import { sameJson, type JsonObject, type JsonValue } from "./json";

export interface RunSettings {
  config: JsonObject;
  inputs: readonly RunInput[];
}

export type SettingDifference =
  | { kind: "setting"; spec: SettingSpec; first: JsonValue; second: JsonValue }
  | { kind: "input"; spec: SettingSpec; first: readonly string[]; second: readonly string[]; sameContent: boolean };

export function settingDifferences(
  specs: readonly SettingSpec[],
  first: RunSettings,
  second: RunSettings,
): SettingDifference[] {
  return specs.flatMap((spec): SettingDifference[] => {
    if (spec.role === "output") {
      return [];
    }

    if (spec.role === "input" || spec.role === "input-template") {
      return inputDifference(spec, first, second);
    }

    const left = settingValue(first.config, spec);
    const right = settingValue(second.config, spec);

    return sameJson(left, right) ? [] : [{ kind: "setting", spec, first: left, second: right }];
  });
}

function inputDifference(spec: SettingSpec, first: RunSettings, second: RunSettings): SettingDifference[] {
  const left = inputsOf(first, spec.key);
  const right = inputsOf(second, spec.key);

  const samePaths = sameJson(
    left.map((input) => input.path),
    right.map((input) => input.path),
  );

  const sameContent = sameJson(
    left.map((input) => input.sha256).toSorted(),
    right.map((input) => input.sha256).toSorted(),
  );

  return samePaths && sameContent
    ? []
    : [
        {
          kind: "input",
          spec,
          first: left.map((input) => input.path),
          second: right.map((input) => input.path),
          sameContent,
        },
      ];
}

function inputsOf(settings: RunSettings, key: string): RunInput[] {
  return settings.inputs.filter((input) => input.setting === key);
}
