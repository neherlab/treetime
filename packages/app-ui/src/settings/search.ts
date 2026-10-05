import type { SparseConfig, SettingSpec } from "@neherlab/app-contracts";

import { fuzzyFilter } from "../fuzzy";
import { isChanged } from "./config";

export function matchingSpecs(
  specs: readonly SettingSpec[],
  config: SparseConfig,
  search: string,
  changedOnly: boolean,
): SettingSpec[] {
  const candidates = changedOnly ? specs.filter((spec) => isChanged(config, spec)) : specs;

  return fuzzyFilter(candidates, search, (spec) => [spec.key, spec.flag, spec.label, spec.help, spec.more].join(" "));
}
