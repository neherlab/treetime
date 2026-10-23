import type { CommandSettings, SettingSpec } from "@neherlab/app-contracts";

export function groupedSpecs(settings: CommandSettings, specs: readonly SettingSpec[]): Array<[string, SettingSpec[]]> {
  return settings.groups.flatMap((group) => {
    const members = specs.filter((spec) => spec.group === group);

    return members.length > 0 ? [[group, members]] : [];
  });
}
