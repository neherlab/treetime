import type { JsonValue, SparseConfig } from "@neherlab/app-contracts";
import type { DeepPartialSkipArrayKey } from "react-hook-form";

import { isJsonObject } from "../settings/json";

type FormScalar = string | number | boolean;

type FormField = FormScalar | FormScalar[];

type FormRecord = { [field: string]: FormField };

type PartialRecord = DeepPartialSkipArrayKey<FormRecord>;

export type FormValue = FormField | FormRecord[];

export type FormConfig = { [key: string]: FormValue | undefined };

export function toFormConfig(config: SparseConfig): FormConfig {
  return Object.fromEntries(Object.entries(config).map(([key, value]) => [key, formValue(value)]));
}

export function toFormValue(value: JsonValue | undefined): FormValue | undefined {
  return value === undefined ? undefined : formValue(value);
}

export function fromFormConfig(values: DeepPartialSkipArrayKey<FormConfig>): SparseConfig {
  return definedEntries(values, (value) => (Array.isArray(value) ? value.map(configItem) : value));
}

function configItem(item: FormScalar | PartialRecord): JsonValue {
  return isPartialRecord(item) ? definedEntries(item, (field) => field) : item;
}

function definedEntries<T>(values: { [key: string]: T | undefined }, convert: (value: T) => JsonValue): SparseConfig {
  return Object.fromEntries(
    Object.entries(values).flatMap(([key, value]) => (value === undefined ? [] : [[key, convert(value)]])),
  );
}

function formValue(value: JsonValue): FormValue {
  if (isFormField(value)) {
    return value;
  }

  if (Array.isArray(value) && value.every(isFormRecord)) {
    return value;
  }

  throw new Error(`The settings form cannot hold the value ${JSON.stringify(value)}`);
}

function isPartialRecord(item: FormScalar | PartialRecord): item is PartialRecord {
  return typeof item === "object";
}

function isFormRecord(value: JsonValue): value is FormRecord {
  return isJsonObject(value) && Object.values(value).every(isFormField);
}

function isFormField(value: JsonValue): value is FormField {
  return isFormScalar(value) || (Array.isArray(value) && value.every(isFormScalar));
}

function isFormScalar(value: JsonValue): value is FormScalar {
  return !Array.isArray(value) && !isJsonObject(value);
}
