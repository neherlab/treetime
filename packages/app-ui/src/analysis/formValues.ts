import type { JsonValue, SparseConfig } from "@neherlab/app-contracts";

import { isJsonObject } from "../settings/json";

type FormScalar = string | number | boolean;

type FormEntry = FormScalar | FormScalar[];

export type FormValue = FormEntry | { [key: string]: FormEntry | undefined };

export type FormConfig = { [key: string]: FormValue | undefined };

export function toFormConfig(config: SparseConfig): FormConfig {
  return Object.fromEntries(Object.entries(config).map(([key, value]) => [key, formValue(value)]));
}

export function toFormValue(value: JsonValue | undefined): FormValue | undefined {
  return value === undefined ? undefined : formValue(value);
}

function formValue(value: JsonValue): FormValue {
  if (isJsonObject(value)) {
    return Object.fromEntries(Object.entries(value).map(([key, child]) => [key, toFormEntry(child)]));
  }

  return toFormEntry(value);
}

export function fromFormConfig(values: FormConfig): SparseConfig {
  return Object.fromEntries(
    Object.entries(values).flatMap(([key, value]) => (value === undefined ? [] : [[key, fromFormValue(value)]])),
  );
}

export function fromFormValue(value: FormValue): JsonValue {
  if (isFormGroup(value)) {
    return Object.fromEntries(
      Object.entries(value).flatMap(([key, child]) => (child === undefined ? [] : [[key, child]])),
    );
  }

  return value;
}

function toFormEntry(value: JsonValue): FormEntry {
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

function isFormGroup(value: FormValue): value is { [key: string]: FormEntry | undefined } {
  return !Array.isArray(value) && typeof value === "object";
}
