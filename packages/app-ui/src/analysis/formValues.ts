import type { JsonValue, SparseConfig } from "@neherlab/app-contracts";

import { isJsonObject } from "../settings/json";

type FormScalar = string | number | boolean;

export type FormValue = FormScalar | FormScalar[];

export type FormConfig = { [key: string]: FormValue | undefined };

export function toFormConfig(config: SparseConfig): FormConfig {
  return Object.fromEntries(Object.entries(config).map(([key, value]) => [key, formValue(value)]));
}

export function toFormValue(value: JsonValue | undefined): FormValue | undefined {
  return value === undefined ? undefined : formValue(value);
}

export function fromFormConfig(values: FormConfig): SparseConfig {
  return Object.fromEntries(
    Object.entries(values).flatMap(([key, value]) => (value === undefined ? [] : [[key, value]])),
  );
}

function formValue(value: JsonValue): FormValue {
  if (isFormScalar(value)) {
    return value;
  }

  if (Array.isArray(value) && value.every(isFormScalar)) {
    return value;
  }

  throw new Error(`The settings form cannot hold the value ${JSON.stringify(value)}`);
}

function isFormScalar(value: JsonValue): value is FormScalar {
  return !Array.isArray(value) && !isJsonObject(value);
}
