import * as z from "zod";

import type { JsonObject, JsonValue } from "../settings/json";

type FormScalar = string | number | boolean | null;

type FormValue = FormScalar | FormScalar[] | Record<string, FormScalar | FormScalar[]>;

export type FormConfig = Record<string, FormValue>;

const zFormScalar = z.union([z.string(), z.number(), z.boolean(), z.null()]);

const zFormValue = z.union([
  zFormScalar,
  z.array(zFormScalar),
  z.record(z.string(), z.union([zFormScalar, z.array(zFormScalar)])),
]);

const zFormConfig = z.record(z.string(), zFormValue);

export function toFormConfig(config: JsonObject): FormConfig {
  return zFormConfig.parse(config);
}

export function toFormValue(value: JsonValue): FormValue {
  return zFormValue.parse(value);
}
