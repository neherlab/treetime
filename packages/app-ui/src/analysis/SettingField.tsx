import type { AppCommand, SettingSpec, SparseConfig } from "@neherlab/app-contracts";
import { Badge } from "@neherlab/design/ui/badge";
import { Button } from "@neherlab/design/ui/button";
import { Label } from "@neherlab/design/ui/label";
import { cn } from "cn";
import { use, useCallback } from "react";
import { useFormContext, useFormState } from "react-hook-form";
import RotateCcw from "~icons/lucide/rotate-ccw";

import { isChanged, resetValue } from "../settings/config";
import { SettingHelpPopover, SettingHelpText } from "../settings/help";
import { settingFieldId } from "./fieldIds";
import { toFormValue, type FormConfig } from "./formValues";
import { SettingControl } from "./SettingControl";
import { WorkspaceContext } from "./workspaceContext";

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
  const commandLine = use(WorkspaceContext)?.check?.command_lines.find((line) => line.key === spec.key);

  const reset = useCallback(
    () => setValue(spec.key, toFormValue(resetValue(spec)), { shouldDirty: true, shouldValidate: true }),
    [setValue, spec],
  );

  return (
    <div
      data-changed={changed || undefined}
      className="grid gap-x-4.5 gap-y-1 border-t py-2.5 pr-3.5 pl-7.5 @xl:grid-cols-[minmax(0,1fr)_16.25rem]"
    >
      <div className="flex flex-wrap items-center gap-1.5">
        <Label htmlFor={settingFieldId(spec.key, "all")}>{label}</Label>
        <SettingHelpPopover spec={spec} commandLine={commandLine} />
        {changed && (
          <Badge variant="secondary" className="text-primary h-4 px-1.5 text-[0.6875rem]">
            changed
          </Badge>
        )}
      </div>
      <div className="flex items-start gap-1.5 @xl:col-start-2 @xl:row-span-2 @xl:row-start-1">
        <div className="grid min-w-0 flex-1 gap-1">
          <SettingControl command={command} spec={spec} label={label} scope="all" />
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
      <SettingHelpText spec={spec} />
    </div>
  );
}
