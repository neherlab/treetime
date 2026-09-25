import type { SettingSpec } from "./catalog";
import { isChanged } from "./config";
import type { JsonObject } from "./json";

export function matchingSpecs(
  specs: readonly SettingSpec[],
  config: JsonObject,
  search: string,
  changedOnly: boolean,
): SettingSpec[] {
  const words = search
    .toLowerCase()
    .split(/\s+/u)
    .filter((word) => word !== "");

  return specs.filter((spec) => {
    if (changedOnly && !isChanged(config, spec)) {
      return false;
    }

    const haystack = [spec.key, spec.flag, spec.label, spec.help, spec.more].join(" ").toLowerCase();

    return words.every((word) => haystack.includes(word));
  });
}
