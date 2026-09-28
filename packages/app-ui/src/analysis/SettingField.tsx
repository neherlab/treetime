import type { AppCommand } from "@neherlab/app-contracts";
import { RotateCcw } from "lucide-react";
import { useCallback } from "react";
import { useFormContext, useFormState } from "react-hook-form";

import type { SettingSpec } from "../settings/catalog";
import { isChanged, resetValue } from "../settings/config";
import { isJsonObject, type JsonObject, type JsonValue } from "../settings/json";
import { formatList } from "../settings/lists";
import { Button } from "../ui/button";
import { cn } from "../ui/cn";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "../ui/collapsible";
import { Label } from "../ui/label";
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
        "grid gap-x-4.5 gap-y-1 border-t py-2.5 pr-3.5 pl-7.5 @xl:grid-cols-[minmax(0,1fr)_16.25rem]",
        changed && "bg-accent/50",
      )}
    >
      <div className="flex flex-wrap items-baseline gap-2">
        <Label htmlFor={settingFieldId(spec.key)}>{label}</Label>
        <code className="text-muted-foreground font-mono text-xs">{spec.flag}</code>
        {spec.role === "setting" && (
          <span className="text-muted-foreground text-xs">Default: {defaultText(spec.default_value)}</span>
        )}
      </div>
      <div className="flex items-start gap-1.5 @xl:col-start-2 @xl:row-span-2 @xl:row-start-1">
        <div className="grid min-w-0 flex-1 gap-1">
          <SettingControl command={command} spec={spec} label={label} />
          {error !== null && <p className="text-destructive text-xs">{error}</p>}
        </div>
        {spec.role === "setting" && (
          <Button
            type="button"
            variant="ghost"
            size="icon-sm"
            onClick={reset}
            title="Reset to default"
            aria-label={`Reset ${label} to default`}
            className={cn(!changed && "invisible")}
          >
            <RotateCcw aria-hidden />
          </Button>
        )}
      </div>
      <SettingHelp spec={spec} />
    </div>
  );
}

export function SettingHelp({ spec }: { spec: SettingSpec }) {
  return (
    <div className="text-muted-foreground max-w-[72ch] text-[0.8125rem]">
      {spec.help}
      {spec.more !== "" && (
        <Collapsible className="inline">
          <CollapsibleTrigger className="text-primary ml-1 underline-offset-4 hover:underline">More</CollapsibleTrigger>
          <CollapsibleContent className="whitespace-pre-line">{spec.more}</CollapsibleContent>
        </Collapsible>
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
