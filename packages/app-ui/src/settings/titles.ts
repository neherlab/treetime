import type { AppCommand, ConfigCode } from "@neherlab/app-contracts";

import type { SettingSpec } from "./catalog";
import { COMMAND_INFO } from "./commands";

const TITLE_FLAGS = 2;

export function autoTitle(command: AppCommand, chips: readonly string[]): string {
  const label = COMMAND_INFO[command].label;

  if (chips.length === 0) {
    return `${label}, defaults`;
  }

  const extra = chips.length > TITLE_FLAGS ? ` +${chips.length - TITLE_FLAGS}` : "";

  return `${label}, ${chips.slice(0, TITLE_FLAGS).join(" ")}${extra}`;
}

export function changedChips(
  code: ConfigCode | null,
  specs: readonly SettingSpec[],
  changed: readonly string[],
): string[] {
  if (code === null) {
    return settingFlags(specs, changed);
  }

  return code.command_line.flatMap((line) => (line.kind === "changed" ? [line.text] : []));
}

export function settingFlags(specs: readonly SettingSpec[], keys: readonly string[]): string[] {
  const flags = new Map(specs.map((spec) => [spec.key, spec.flag]));

  return keys.map((key) => flags.get(key) ?? key);
}
