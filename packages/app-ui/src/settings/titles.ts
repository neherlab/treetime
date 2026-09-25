import type { AppCommand } from "@neherlab/app-contracts";

import { flagTokens } from "./commandLine";
import { COMMAND_INFO } from "./commands";
import { changedSpecs, settingValue } from "./config";
import type { JsonObject } from "./json";
import type { SettingSpec } from "./schema";

const TITLE_FLAGS = 2;

export function autoTitle(command: AppCommand, specs: readonly SettingSpec[], config: JsonObject): string {
  const chips = changedSpecs(specs, config).map((spec) => {
    const tokens = flagTokens(spec, settingValue(config, spec));

    return tokens === null ? spec.flag : tokens.join(" ");
  });

  const label = COMMAND_INFO[command].label;

  if (chips.length === 0) {
    return `${label}, defaults`;
  }

  const extra = chips.length > TITLE_FLAGS ? ` +${chips.length - TITLE_FLAGS}` : "";

  return `${label}, ${chips.slice(0, TITLE_FLAGS).join(" ")}${extra}`;
}

export function settingFlags(specs: readonly SettingSpec[], keys: readonly string[]): string[] {
  const flags = new Map(specs.map((spec) => [spec.key, spec.flag]));

  return keys.map((key) => flags.get(key) ?? key);
}
