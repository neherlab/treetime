export function settingFieldId(key: string): string {
  return `setting-${key.replaceAll(".", "-")}`;
}
