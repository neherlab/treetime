import type { RunRecord, SparseConfig } from "@neherlab/app-contracts";

import type { Locale } from "../i18n/locale";
import { commandSettings } from "../settings/catalog";
import { effectiveValue } from "../settings/config";
import { defaultText } from "../settings/help";

export function settingText(record: RunRecord, key: string, locale: Locale): string {
  const config: SparseConfig = record.config;
  const spec = commandSettings(record.command).settings.find((candidate) => candidate.key === key);

  return spec === undefined
    ? defaultText(config[key], locale)
    : defaultText(effectiveValue(config, spec), locale, spec.options);
}
