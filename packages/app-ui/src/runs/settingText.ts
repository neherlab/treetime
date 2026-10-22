import type { RunRecord, SparseConfig } from "@neherlab/app-contracts";

import { defaultText } from "../analysis/SettingField";
import { commandSettings } from "../settings/catalog";
import { effectiveValue } from "../settings/config";

export function settingText(record: RunRecord, key: string): string {
  const config: SparseConfig = record.config;
  const spec = commandSettings(record.command).settings.find((candidate) => candidate.key === key);

  return defaultText(spec === undefined ? config[key] : effectiveValue(config, spec));
}
