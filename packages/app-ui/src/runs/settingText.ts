import type { RunRecord, SparseConfig } from "@neherlab/app-contracts";

import { defaultText } from "../analysis/SettingField";

export function settingText(record: RunRecord, key: string): string {
  const config: SparseConfig = record.config;

  return defaultText(config[key]);
}
