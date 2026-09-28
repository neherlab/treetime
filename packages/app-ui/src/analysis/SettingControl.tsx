import type { AppCommand, SettingKind, SettingOption } from "@neherlab/app-contracts";
import { useCallback, useMemo, useState } from "react";
import { useController } from "react-hook-form";

import type { SettingSpec } from "../settings/catalog";
import { baseName, pathList } from "../settings/inputs";
import { isJsonObject, sameJson, type JsonValue } from "../settings/json";
import { formatList, parseList } from "../settings/lists";
import { cn } from "../ui/cn";
import { Input } from "../ui/input";
import { NativeSelect, NativeSelectOption } from "../ui/native-select";
import { Switch } from "../ui/switch";
import { ToggleGroup, ToggleGroupItem } from "../ui/toggle-group";
import { settingFieldId } from "./fieldIds";
import { toFormValue, type FormConfig } from "./formValues";
import { PathPicker } from "./PathPicker";

const CONTROL_HEIGHT = "h-8";

const NO_EXTENSIONS: readonly string[] = [];

interface ControlProps {
  command: AppCommand;
  spec: SettingSpec;
  label: string;
  className: string | undefined;
}

const CONTROLS: Record<SettingKind, (props: ControlProps) => React.ReactNode> = {
  switch: SwitchControl,
  tristate: TristateControl,
  enum: EnumControl,
  integer: NumberControl,
  number: NumberControl,
  text: TextControl,
  list: ListControl,
  "enum-list": EnumListControl,
};

export function SettingControl({
  command,
  spec,
  label,
  className,
}: {
  command: AppCommand;
  spec: SettingSpec;
  label: string;
  className?: string | undefined;
}) {
  if (spec.role === "output") {
    return (
      <Input
        id={settingFieldId(spec.key)}
        type="text"
        disabled
        aria-label={label}
        placeholder="Set by the app for each run"
        className={cn(CONTROL_HEIGHT, className)}
      />
    );
  }

  if (spec.role === "input") {
    return <PathControl command={command} spec={spec} label={label} className={className} />;
  }

  const Control = CONTROLS[spec.kind];

  return <Control command={command} spec={spec} label={label} className={className} />;
}

function PathControl({ command, spec, label, className }: ControlProps) {
  const { value } = useSetting(spec);
  const paths = pathList(value);

  return (
    <div className={cn("flex min-w-0 items-center gap-2", className)}>
      <span
        id={settingFieldId(spec.key)}
        className="min-w-0 flex-1 truncate font-mono text-xs"
        title={paths.join("\n")}
      >
        {paths.length > 0 ? paths.map(baseName).join(", ") : "No file"}
      </span>
      <PathPicker
        command={command}
        settingKey={spec.key}
        title={label}
        extensions={NO_EXTENSIONS}
        list={spec.kind === "list"}
        filled={paths.length > 0}
      />
    </div>
  );
}

function SwitchControl({ spec, label, className }: ControlProps) {
  const { value, set } = useSetting(spec);

  return (
    <Switch
      id={settingFieldId(spec.key)}
      aria-label={label}
      checked={value === true}
      onCheckedChange={set}
      className={className}
    />
  );
}

function TristateControl({ spec, label, className }: ControlProps) {
  const { value, set } = useSetting(spec);

  const onChange = useCallback(
    (event: React.ChangeEvent<HTMLSelectElement>) =>
      set(event.target.value === "" ? null : event.target.value === "true"),
    [set],
  );

  return (
    <NativeSelect
      id={settingFieldId(spec.key)}
      size="sm"
      aria-label={label}
      value={value === null ? "" : scalarText(value)}
      onChange={onChange}
      className={className}
    >
      <NativeSelectOption value="">Automatic</NativeSelectOption>
      <NativeSelectOption value="true">On</NativeSelectOption>
      <NativeSelectOption value="false">Off</NativeSelectOption>
    </NativeSelect>
  );
}

function EnumControl({ spec, label, className }: ControlProps) {
  const { value, set, invalid } = useSetting(spec);

  const onChange = useCallback(
    (event: React.ChangeEvent<HTMLSelectElement>) => set(event.target.value === "" ? null : event.target.value),
    [set],
  );

  return (
    <NativeSelect
      id={settingFieldId(spec.key)}
      size="sm"
      aria-label={label}
      aria-invalid={invalid}
      value={value === null ? "" : scalarText(value)}
      onChange={onChange}
      className={className}
    >
      {spec.nullable && <NativeSelectOption value="">Not set</NativeSelectOption>}
      {spec.options.map((option) => (
        <NativeSelectOption key={option.value} value={option.value} title={option.help}>
          {option.value}
        </NativeSelectOption>
      ))}
    </NativeSelect>
  );
}

function NumberControl({ spec, label, className }: ControlProps) {
  const { value, set, onBlur, invalid } = useSetting(spec);

  const onChange = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => set(event.target.value === "" ? null : Number(event.target.value)),
    [set],
  );

  return (
    <Input
      id={settingFieldId(spec.key)}
      type="number"
      aria-label={label}
      aria-invalid={invalid}
      step={spec.kind === "integer" ? 1 : "any"}
      min={spec.minimum ?? undefined}
      value={value === null ? "" : scalarText(value)}
      placeholder={spec.nullable ? "Not set" : ""}
      onChange={onChange}
      onBlur={onBlur}
      className={cn(CONTROL_HEIGHT, "font-mono", className)}
    />
  );
}

function TextControl({ spec, label, className }: ControlProps) {
  const { value, set, onBlur, invalid } = useSetting(spec);

  const onChange = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) =>
      set(event.target.value === "" && spec.nullable ? null : event.target.value),
    [set, spec.nullable],
  );

  return (
    <Input
      id={settingFieldId(spec.key)}
      type="text"
      aria-label={label}
      aria-invalid={invalid}
      value={value === null ? "" : scalarText(value)}
      placeholder={spec.nullable ? "Not set" : ""}
      onChange={onChange}
      onBlur={onBlur}
      className={cn(CONTROL_HEIGHT, className)}
    />
  );
}

function ListControl({ spec, label, className }: ControlProps) {
  const { value, set, invalid } = useSetting(spec);
  const [text, setText] = useState(() => formatList(value));
  const shown = sameJson(parseList(text, spec.item_kind), value) ? text : formatList(value);

  const onChange = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => {
      setText(event.target.value);
      set(parseList(event.target.value, spec.item_kind));
    },
    [set, spec.item_kind],
  );

  return (
    <Input
      id={settingFieldId(spec.key)}
      type="text"
      aria-label={label}
      aria-invalid={invalid}
      value={shown}
      placeholder={spec.item_kind === "string" ? "Values, separated by spaces" : "Numbers, separated by spaces"}
      onChange={onChange}
      className={cn(CONTROL_HEIGHT, "font-mono", className)}
    />
  );
}

function EnumListControl({ spec, label, className }: ControlProps) {
  const { value, set } = useSetting(spec);
  const selected = useMemo(() => (Array.isArray(value) ? value.map(scalarText) : []), [value]);

  const onValueChange = useCallback(
    (next: string[]) => set(spec.options.flatMap((option) => (next.includes(option.value) ? [option.value] : []))),
    [set, spec.options],
  );

  return (
    <ToggleGroup
      id={settingFieldId(spec.key)}
      aria-label={label}
      multiple
      variant="outline"
      size="sm"
      spacing={1}
      value={selected}
      onValueChange={onValueChange}
      className={cn("flex-wrap", className)}
    >
      {spec.options.map((option) => (
        <EnumListOption key={option.value} option={option} />
      ))}
    </ToggleGroup>
  );
}

function EnumListOption({ option }: { option: SettingOption }) {
  return (
    <ToggleGroupItem value={option.value} title={option.help} className="h-7 px-2 font-mono text-xs">
      {option.value}
    </ToggleGroupItem>
  );
}

function scalarText(value: JsonValue): string {
  return Array.isArray(value) || isJsonObject(value) ? JSON.stringify(value) : `${value}`;
}

function useSetting(spec: SettingSpec) {
  const { field, fieldState } = useController<FormConfig>({ name: spec.key });
  const value: JsonValue = field.value ?? null;
  const { onChange } = field;
  const set = useCallback((next: JsonValue) => onChange(toFormValue(next)), [onChange]);

  return { value, set, onBlur: field.onBlur, invalid: fieldState.error !== undefined };
}
