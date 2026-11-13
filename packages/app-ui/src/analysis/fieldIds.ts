export function settingFieldId(key: string, scope: FieldScope = "all"): string {
  return `setting-${scope}-${key.replaceAll(".", "-")}`;
}

export type FieldScope = "all" | "main" | "change";
