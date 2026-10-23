import type { AppCommand, SettingSpec, SparseConfig } from "@neherlab/app-contracts";
import { cn } from "cn";
import { useCallback } from "react";
import { useFormContext, useFormState } from "react-hook-form";
import RotateCcw from "~icons/lucide/rotate-ccw";

import { isChanged, resetValue } from "../settings/config";
import { SettingHelp, SettingSummary } from "../settings/help";
import { Button } from "../ui/button";
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
  config: SparseConfig;
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
      <SettingSummary spec={spec} label={<Label htmlFor={settingFieldId(spec.key)}>{label}</Label>} />
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
