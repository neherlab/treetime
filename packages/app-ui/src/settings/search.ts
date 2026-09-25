import { wordMatcher } from "../text";
import type { SettingSpec } from "./catalog";
import { isChanged } from "./config";
import type { JsonObject } from "./json";

export function matchingSpecs(
  specs: readonly SettingSpec[],
  config: JsonObject,
  search: string,
  changedOnly: boolean,
): SettingSpec[] {
  const matches = wordMatcher(search);

  return specs.filter(
    (spec) =>
      (!changedOnly || isChanged(config, spec)) &&
      matches([spec.key, spec.flag, spec.label, spec.help, spec.more].join(" ")),
  );
}
