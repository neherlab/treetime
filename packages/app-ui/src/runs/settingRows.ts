import type { SettingSpec, SparseConfig } from "@neherlab/app-contracts";

export type ShownSpec = Pick<SettingSpec, "key" | "role" | "set_by_app">;

export function shownSettings<T extends ShownSpec>(
  specs: readonly T[],
  config: SparseConfig,
  changed: ReadonlySet<string>,
  changedOnly: boolean,
): T[] {
  return specs.filter(
    (spec) =>
      !spec.set_by_app &&
      (!changedOnly || changed.has(spec.key) || (spec.role === "input" && config[spec.key] !== undefined)),
  );
}
