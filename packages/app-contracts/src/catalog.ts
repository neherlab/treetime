import catalogJson from "../setting-catalog.json";
import type { SettingCatalog } from "./generated/types.gen";
import { zSettingCatalog } from "./generated/zod.gen";

const catalog: unknown = catalogJson;

assertSettingCatalog(catalog);

export const settingCatalog: SettingCatalog = catalog;

function assertSettingCatalog(value: unknown): asserts value is SettingCatalog {
  zSettingCatalog.parse(value);
}
