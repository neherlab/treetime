import { fuzzyFilter } from "../fuzzy";
import type { SettingSpec } from "./catalog";
import { isChanged } from "./config";
import type { JsonObject } from "./json";

export function matchingSpecs(
  specs: readonly SettingSpec[],
  config: JsonObject,
  search: string,
  changedOnly: boolean,
): SettingSpec[] {
  const candidates = changedOnly ? specs.filter((spec) => isChanged(config, spec)) : specs;

  return fuzzyFilter(candidates, search, (spec) => [spec.key, spec.flag, spec.label, spec.help, spec.more].join(" "));
}
