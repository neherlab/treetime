import type { CommandSettings, SettingPurpose, SettingSpec } from "@neherlab/app-contracts";

export function groupedSpecs(settings: CommandSettings, specs: readonly SettingSpec[]): Array<[string, SettingSpec[]]> {
  return settings.groups.flatMap((group) => {
    const members = specs.filter((spec) => spec.group === group);

    return members.length > 0 ? [[group, members]] : [];
  });
}

export function settingSections(settings: CommandSettings, specs: readonly SettingSpec[]): SettingSectionView[] {
  const offered = specs.filter((spec) => !spec.set_by_app);

  return settings.sections.flatMap(({ purpose, title, groups }) => {
    const members = groups.flatMap((group) => {
      const inGroup = offered.filter((spec) => spec.purpose === purpose && spec.group === group);

      return inGroup.length > 0 ? [{ group, specs: inGroup }] : [];
    });

    return members.length > 0 ? [{ purpose, title, groups: members }] : [];
  });
}

export interface SettingSectionView {
  purpose: SettingPurpose;
  title: string;
  groups: SettingGroupView[];
}

export interface SettingGroupView {
  group: string;
  specs: SettingSpec[];
}
