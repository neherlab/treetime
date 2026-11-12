import type { SettingSpec, SparseConfig } from "@neherlab/app-contracts";

import { normalizeConfig } from "./config";
import { withSetting } from "./json";

export function withRunSettings(
  specs: readonly SettingSpec[],
  config: SparseConfig,
  changedKeys: readonly string[],
  runSettings: SparseConfig,
): SparseConfig {
  const changed = new Set(changedKeys);

  const applied = specs
    .filter((spec) => spec.role !== "input" && !spec.set_by_app && changed.has(spec.key))
    .reduce((next, spec) => withSetting(next, spec.key, runSettings[spec.key]), config);

  return normalizeConfig(specs, applied);
}
