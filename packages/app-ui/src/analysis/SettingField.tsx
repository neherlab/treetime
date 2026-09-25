import type { AppCommand } from "@neherlab/app-contracts";
import { useCallback } from "react";
import { useFormContext, useFormState } from "react-hook-form";

import type { SettingSpec } from "../settings/catalog";
import { isChanged, resetValue } from "../settings/config";
import { isJsonObject, type JsonObject, type JsonValue } from "../settings/json";
import { formatList } from "../settings/lists";
import { cn } from "../ui";
import { settingFieldId } from "./fieldIds";
import { toFormValue, type FormConfig } from "./formValues";
import { SettingControl } from "./SettingControl";

export function SettingField({
  command,
  spec,
  config,
}: {
  command: AppCommand;
  spec: SettingSpec;
  config: JsonObject;
}) {
  const { setValue, getFieldState } = useFormContext<FormConfig>();
  const formState = useFormState<FormConfig>({ name: spec.key });
  const changed = isChanged(config, spec);
  const label = spec.label;
  const error = getFieldState(spec.key, formState).error?.message ?? null;

  const reset = useCallback(
    () => setValue(spec.key, toFormValue(resetValue(spec)), { shouldDirty: true, shouldValidate: true }),
    [setValue, spec],
  );

  return (
    <div
      className={cn(
        "border-surface-3 grid gap-x-4.5 gap-y-1 border-t py-2 pr-3.5 pl-7.5 md:grid-cols-[minmax(0,1fr)_16.25rem]",
        changed && "from-accent-subtle bg-gradient-to-r to-transparent to-60%",
      )}
    >
      <div className="flex flex-wrap items-baseline gap-2">
        <label htmlFor={settingFieldId(spec.key)} className="font-bold">
          {label}
        </label>
        <code className="text-ink-faint font-mono text-xs">{spec.flag}</code>
        {spec.role === "setting" && (
          <span className="text-ink-faint text-xs">Default: {defaultText(spec.default_value)}</span>
        )}
      </div>
      <div className="flex items-start gap-1.5 md:col-start-2 md:row-span-2 md:row-start-1">
        <div className="min-w-0 flex-1">
          <SettingControl command={command} spec={spec} label={label} />
          {error !== null && <p className="text-signal-danger mt-0.5 text-xs">{error}</p>}
        </div>
        {spec.role === "setting" && (
          <button
            type="button"
            onClick={reset}
            title="Reset to default"
            className={cn(
              "text-ink-muted hover:bg-surface-2 rounded-md px-2 py-1 text-xs font-bold",
              !changed && "invisible",
            )}
          >
            Reset
          </button>
        )}
      </div>
      <SettingHelp spec={spec} />
    </div>
  );
}

export function SettingHelp({ spec }: { spec: SettingSpec }) {
  return (
    <div className="text-ink-muted max-w-[72ch] text-[0.8125rem]">
      {spec.help}
      {spec.more !== "" && (
        <details className="inline">
          <summary className="text-accent ml-1 inline cursor-pointer">More</summary>
          <span className="block whitespace-pre-line">{spec.more}</span>
        </details>
      )}
    </div>
  );
}

export function defaultText(value: JsonValue): string {
  if (value === null) {
    return "not set";
  }

  if (Array.isArray(value)) {
    return value.length === 0 ? "none" : formatList(value);
  }

  if (value === true) {
    return "on";
  }

  if (value === false) {
    return "off";
  }

  return isJsonObject(value) ? JSON.stringify(value) : String(value);
}
