import type { AppCommand, SettingKind, SettingOption } from "@neherlab/app-contracts";
import { useCallback, useMemo, useState } from "react";
import { useController } from "react-hook-form";

import type { SettingSpec } from "../settings/catalog";
import { baseName, pathList } from "../settings/inputs";
import { isJsonObject, sameJson, type JsonValue } from "../settings/json";
import { formatList, parseList } from "../settings/lists";
import { Switch, cn } from "../ui";
import { settingFieldId } from "./fieldIds";
import { toFormValue, type FormConfig } from "./formValues";
import { PathPicker } from "./PathPicker";

const INPUT_CLASS =
  "border-line-strong bg-surface-1 focus-visible:ring-accent aria-[invalid=true]:border-signal-danger w-full rounded-md border px-2 py-1 outline-none focus-visible:ring-2 disabled:opacity-60";

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
      <input
        id={settingFieldId(spec.key)}
        type="text"
        disabled
        aria-label={label}
        placeholder="Set by the app for each run"
        className={cn(INPUT_CLASS, className)}
      />
    );
  }

  if (spec.role !== "setting") {
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
      label={label}
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
    <select
      id={settingFieldId(spec.key)}
      aria-label={label}
      value={value === null ? "" : scalarText(value)}
      onChange={onChange}
      className={cn(INPUT_CLASS, className)}
    >
      <option value="">Automatic</option>
      <option value="true">On</option>
      <option value="false">Off</option>
    </select>
  );
}

function EnumControl({ spec, label, className }: ControlProps) {
  const { value, set, invalid } = useSetting(spec);

  const onChange = useCallback(
    (event: React.ChangeEvent<HTMLSelectElement>) => set(event.target.value === "" ? null : event.target.value),
    [set],
  );

  return (
    <select
      id={settingFieldId(spec.key)}
      aria-label={label}
      aria-invalid={invalid}
      value={value === null ? "" : scalarText(value)}
      onChange={onChange}
      className={cn(INPUT_CLASS, className)}
    >
      {spec.nullable && <option value="">Not set</option>}
      {spec.options.map((option) => (
        <option key={option.value} value={option.value} title={option.help}>
          {option.value}
        </option>
      ))}
    </select>
  );
}

function NumberControl({ spec, label, className }: ControlProps) {
  const { value, set, onBlur, invalid } = useSetting(spec);

  const onChange = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => set(event.target.value === "" ? null : Number(event.target.value)),
    [set],
  );

  return (
    <input
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
      className={cn(INPUT_CLASS, className)}
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
    <input
      id={settingFieldId(spec.key)}
      type="text"
      aria-label={label}
      aria-invalid={invalid}
      value={value === null ? "" : scalarText(value)}
      placeholder={spec.nullable ? "Not set" : ""}
      onChange={onChange}
      onBlur={onBlur}
      className={cn(INPUT_CLASS, className)}
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
    <input
      id={settingFieldId(spec.key)}
      type="text"
      aria-label={label}
      aria-invalid={invalid}
      value={shown}
      placeholder={spec.item_kind === "string" ? "Values, separated by spaces" : "Numbers, separated by spaces"}
      onChange={onChange}
      className={cn(INPUT_CLASS, className)}
    />
  );
}

function EnumListControl({ spec, label, className }: ControlProps) {
  const { value, set } = useSetting(spec);
  const selected = useMemo(() => (Array.isArray(value) ? value : []), [value]);

  return (
    <fieldset id={settingFieldId(spec.key)} className={cn("m-0 flex flex-wrap gap-1 border-0 p-0", className)}>
      <legend className="sr-only">{label}</legend>
      {spec.options.map((option) => (
        <EnumListOption key={option.value} spec={spec} option={option} selected={selected} set={set} />
      ))}
    </fieldset>
  );
}

function EnumListOption({
  spec,
  option,
  selected,
  set,
}: {
  spec: SettingSpec;
  option: SettingOption;
  selected: readonly JsonValue[];
  set: (value: JsonValue) => void;
}) {
  const onChange = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => {
      const chosen = new Set(selected);

      if (event.target.checked) {
        chosen.add(option.value);
      } else {
        chosen.delete(option.value);
      }

      set(spec.options.flatMap((candidate) => (chosen.has(candidate.value) ? [candidate.value] : [])));
    },
    [option.value, selected, set, spec.options],
  );

  return (
    <label title={option.help} className="border-line inline-flex items-center gap-1 rounded-sm border px-1.5 text-xs">
      <input type="checkbox" checked={selected.includes(option.value)} onChange={onChange} />
      {option.value}
    </label>
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
