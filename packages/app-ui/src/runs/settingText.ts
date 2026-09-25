import type { RunRecordResult } from "@neherlab/app-contracts";

import { defaultText } from "../analysis/SettingField";
import { zJsonObject } from "../settings/json";

export function settingText(record: RunRecordResult, key: string): string {
  return defaultText(zJsonObject.parse(record.config)[key] ?? null);
}
