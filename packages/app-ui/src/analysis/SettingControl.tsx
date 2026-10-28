import type {
  AppCommand,
  JsonValue,
  RecordFieldSpec,
  SettingKind,
  SettingOption,
  SettingSpec,
  SparseConfig,
} from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { Input } from "@neherlab/design/ui/input";
import { NativeSelect, NativeSelectOption } from "@neherlab/design/ui/native-select";
import { Switch } from "@neherlab/design/ui/switch";
import { ToggleGroup, ToggleGroupItem } from "@neherlab/design/ui/toggle-group";
import { Tooltip, TooltipContent, TooltipTrigger } from "@neherlab/design/ui/tooltip";
import { cn } from "cn";
import { useCallback, useMemo, useState } from "react";
import { GridList, GridListItem } from "react-aria-components";
import { useController, useWatch } from "react-hook-form";
import Plus from "~icons/lucide/plus";
import Trash from "~icons/lucide/trash-2";

import { baseName, pathList } from "../settings/inputs";
import { isJsonObject, sameJson, withSetting } from "../settings/json";
import { formatList, parseList } from "../settings/lists";
import { settingFieldId } from "./fieldIds";
import { toFormValue, type FormConfig } from "./formValues";
import { NumberInput } from "./NumberInput";
import { PathPicker } from "./PathPicker";

const CONTROL_HEIGHT = "h-8";

const NO_EXTENSIONS: readonly string[] = [];

const FREE_VALUE = "#value";

const NO_FIELDS: RecordFieldSpec[] = [];

interface ControlProps {
  command: AppCommand;
  spec: SettingSpec;
  label: string;
  className: string | undefined;
}

type EditedSetting = Pick<
  SettingSpec,
  "nullable" | "options" | "item_kind" | "examples" | "documented_default" | "fields"
>;

interface EditorProps {
  id: string | undefined;
  label: string;
  setting: EditedSetting;
  value: JsonValue | undefined;
  set: (next: JsonValue) => void;
  unset: () => void;
  onKeyDown: ((event: React.KeyboardEvent) => void) | undefined;
  onBlur: (() => void) | undefined;
  invalid: boolean;
  className: string | undefined;
}

const EDITORS: Record<SettingKind, (props: EditorProps) => React.ReactNode> = {
  switch: SwitchEditor,
  enum: EnumEditor,
  "enum-or-text": EnumOrTextEditor,
  integer: NumberEditor,
  number: NumberEditor,
  text: TextEditor,
  list: ListEditor,
  "enum-list": EnumListEditor,
  records: RecordsEditor,
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

  return <ValueControl command={command} spec={spec} label={label} className={className} />;
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

function ValueControl({ spec, label, className }: ControlProps) {
  const { value, set, unset, onBlur, invalid } = useSetting(spec);
  const Editor = EDITORS[spec.kind];

  return (
    <Editor
      id={settingFieldId(spec.key)}
      label={label}
      setting={spec}
      value={value}
      set={set}
      unset={unset}
      onBlur={onBlur}
      onKeyDown={undefined}
      invalid={invalid}
      className={className}
    />
  );
}

function SwitchEditor({ id, label, value, set, onKeyDown, className }: EditorProps) {
  return (
    <Switch
      id={id}
      aria-label={label}
      checked={value === true}
      onCheckedChange={set}
      onKeyDown={onKeyDown}
      className={className}
    />
  );
}

function EnumEditor({ id, label, setting, value, set, unset, invalid, onKeyDown, className }: EditorProps) {
  const onChange = useCallback(
    (event: React.ChangeEvent<HTMLSelectElement>) => (event.target.value === "" ? unset() : set(event.target.value)),
    [set, unset],
  );

  return (
    <NativeSelect
      id={id}
      size="sm"
      aria-label={label}
      aria-invalid={invalid}
      value={value === undefined ? "" : scalarText(value)}
      onChange={onChange}
      onKeyDown={onKeyDown}
      className={className}
    >
      {setting.nullable && <NativeSelectOption value="">{unsetText(setting)}</NativeSelectOption>}
      {setting.options.map((option) => (
        <NativeSelectOption key={option.value} value={option.value} title={option.help}>
          {option.value}
        </NativeSelectOption>
      ))}
    </NativeSelect>
  );
}

function EnumOrTextEditor({
  id,
  label,
  setting,
  value,
  set,
  unset,
  onBlur,
  invalid,
  onKeyDown,
  className,
}: EditorProps) {
  const text = value === undefined ? "" : scalarText(value);
  const isOption = setting.options.some((option) => option.value === text);
  const selected = value === undefined ? "" : isOption ? text : FREE_VALUE;
  const example = setting.examples[0];

  const onSelect = useCallback(
    (event: React.ChangeEvent<HTMLSelectElement>) => {
      const next = event.target.value;

      if (next === "") {
        unset();
      } else if (next === FREE_VALUE) {
        set(example ?? "");
      } else {
        set(next);
      }
    },
    [example, set, unset],
  );

  const onText = useCallback((event: React.ChangeEvent<HTMLInputElement>) => set(event.target.value), [set]);

  return (
    <div className={cn("flex min-w-0 items-center gap-2", className)}>
      <NativeSelect
        id={id}
        size="sm"
        aria-label={label}
        aria-invalid={invalid}
        value={selected}
        onChange={onSelect}
        onKeyDown={onKeyDown}
      >
        {setting.nullable && <NativeSelectOption value="">{unsetText(setting)}</NativeSelectOption>}
        {setting.options.map((option) => (
          <NativeSelectOption key={option.value} value={option.value} title={option.help}>
            {option.value}
          </NativeSelectOption>
        ))}
        <NativeSelectOption value={FREE_VALUE}>Value</NativeSelectOption>
      </NativeSelect>
      {selected === FREE_VALUE && (
        <Input
          type="text"
          aria-label={`${label} value`}
          aria-invalid={invalid}
          value={text}
          placeholder={example === undefined ? "" : `For example ${scalarText(example)}`}
          onChange={onText}
          onKeyDown={onKeyDown}
          onBlur={onBlur}
          className={cn(CONTROL_HEIGHT, "w-28 font-mono")}
        />
      )}
    </div>
  );
}

function NumberEditor({ id, label, setting, value, set, unset, onBlur, invalid, onKeyDown, className }: EditorProps) {
  const onValueChange = useCallback(
    (next: number | string | undefined) => (next === undefined ? unset() : set(next)),
    [set, unset],
  );

  return (
    <NumberInput
      id={id}
      aria-label={label}
      aria-invalid={invalid}
      value={value}
      placeholder={setting.nullable ? unsetText(setting) : ""}
      onValueChange={onValueChange}
      onKeyDown={onKeyDown}
      onBlur={onBlur}
      className={cn(CONTROL_HEIGHT, "font-mono", className)}
    />
  );
}

function TextEditor({ id, label, setting, value, set, unset, onBlur, invalid, onKeyDown, className }: EditorProps) {
  const onChange = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) =>
      event.target.value === "" && setting.nullable ? unset() : set(event.target.value),
    [set, setting.nullable, unset],
  );

  return (
    <Input
      id={id}
      type="text"
      aria-label={label}
      aria-invalid={invalid}
      value={value === undefined ? "" : scalarText(value)}
      placeholder={setting.nullable ? unsetText(setting) : ""}
      onChange={onChange}
      onKeyDown={onKeyDown}
      onBlur={onBlur}
      className={cn(CONTROL_HEIGHT, className)}
    />
  );
}

function ListEditor({ id, label, setting, value, set, unset, invalid, onKeyDown, className }: EditorProps) {
  const [text, setText] = useState(() => formatList(value));
  const shown = sameJson(parseList(text, setting.item_kind), value) ? text : formatList(value);

  const onChange = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => {
      setText(event.target.value);
      const items = parseList(event.target.value, setting.item_kind);

      if (items.length === 0 && setting.nullable) {
        unset();
      } else {
        set(items);
      }
    },
    [set, setting.item_kind, setting.nullable, unset],
  );

  return (
    <Input
      id={id}
      type="text"
      aria-label={label}
      aria-invalid={invalid}
      value={shown}
      placeholder={setting.item_kind === "string" ? "Values, separated by spaces" : "Numbers, separated by spaces"}
      onChange={onChange}
      onKeyDown={onKeyDown}
      className={cn(CONTROL_HEIGHT, "font-mono", className)}
    />
  );
}

function EnumListEditor({ id, label, setting, value, set, onKeyDown, className }: EditorProps) {
  const selected = useMemo(() => (Array.isArray(value) ? value.map(scalarText) : []), [value]);

  const onValueChange = useCallback(
    (next: string[]) => set(setting.options.flatMap((option) => (next.includes(option.value) ? [option.value] : []))),
    [set, setting.options],
  );

  return (
    <ToggleGroup
      id={id}
      aria-label={label}
      multiple
      variant="outline"
      size="sm"
      spacing={1}
      value={selected}
      onValueChange={onValueChange}
      onKeyDown={onKeyDown}
      className={cn("flex-wrap", className)}
    >
      {setting.options.map((option) => (
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

function unsetText(setting: EditedSetting): string {
  return setting.documented_default === undefined ? "Not set" : `Default: ${scalarText(setting.documented_default)}`;
}

function RecordsEditor({ id, label, setting, value, set, unset, invalid, className }: EditorProps) {
  const records = useMemo(() => (Array.isArray(value) ? value.filter(isJsonObject) : []), [value]);
  const rows = useMemo(() => records.map((record, index) => ({ id: index, record })), [records]);

  const update = useCallback(
    (next: SparseConfig[]) => (next.length === 0 && setting.nullable ? unset() : set(next)),
    [set, setting.nullable, unset],
  );

  const setField = useCallback(
    (index: number, key: string, next: JsonValue | undefined) =>
      update(records.map((record, position) => (position === index ? withSetting(record, key, next) : record))),
    [records, update],
  );

  const remove = useCallback(
    (index: number) => update(records.filter((_, position) => position !== index)),
    [records, update],
  );

  const add = useCallback(() => update([...records, {}]), [records, update]);

  const dependencies = useMemo(
    () => [setting.fields, label, setField, remove],
    [setting.fields, label, setField, remove],
  );

  const renderRow = useCallback(
    (row: RecordRow) => (
      <RecordItem row={row} fields={setting.fields} label={label} setField={setField} remove={remove} />
    ),
    [label, remove, setField, setting.fields],
  );

  return (
    <div className={cn("flex min-w-0 flex-col gap-2", className)}>
      <GridList
        {...(id === undefined ? {} : { id })}
        aria-label={label}
        items={rows}
        dependencies={dependencies}
        layout="stack"
        keyboardNavigationBehavior="tab"
        disallowTypeAhead
        renderEmptyState={renderNoEntries}
        data-invalid={invalid || undefined}
        className="data-invalid:border-destructive flex flex-col gap-2 rounded-md border p-2 outline-none"
      >
        {renderRow}
      </GridList>
      <Button type="button" variant="outline" size="sm" onClick={add} className="self-start">
        <Plus aria-hidden />
        Add entry
      </Button>
    </div>
  );
}

function RecordItem({
  row,
  fields,
  label,
  setField,
  remove,
}: {
  row: RecordRow;
  fields: RecordFieldSpec[];
  label: string;
  setField: (index: number, key: string, next: JsonValue | undefined) => void;
  remove: (index: number) => void;
}) {
  const entry = `${label} entry ${row.id + 1}`;
  const onRemove = useCallback(() => remove(row.id), [remove, row.id]);

  return (
    <GridListItem
      id={row.id}
      textValue={entry}
      className="focus-visible:ring-ring/50 flex items-end gap-2 rounded-md outline-none focus-visible:ring-3"
    >
      {fields.map((field) => (
        <RecordField key={field.key} field={field} entry={entry} row={row} setField={setField} />
      ))}
      <Tooltip>
        <TooltipTrigger
          render={
            <Button type="button" variant="ghost" size="icon-sm" aria-label={`Remove ${entry}`} onClick={onRemove}>
              <Trash aria-hidden />
            </Button>
          }
        />
        <TooltipContent>Remove this entry</TooltipContent>
      </Tooltip>
    </GridListItem>
  );
}

function RecordField({
  field,
  entry,
  row,
  setField,
}: {
  field: RecordFieldSpec;
  entry: string;
  row: RecordRow;
  setField: (index: number, key: string, next: JsonValue | undefined) => void;
}) {
  const setting = useMemo(() => ({ ...field, fields: NO_FIELDS }), [field]);
  const set = useCallback((next: JsonValue) => setField(row.id, field.key, next), [field.key, row.id, setField]);
  const unset = useCallback(() => setField(row.id, field.key, undefined), [field.key, row.id, setField]);
  const Editor = EDITORS[field.kind];

  return (
    <div className="flex min-w-0 flex-1 flex-col gap-1">
      <span className="text-muted-foreground text-xs" title={field.help}>
        {field.label}
      </span>
      <Editor
        id={undefined}
        label={`${field.label} of ${entry}`}
        setting={setting}
        value={row.record[field.key]}
        set={set}
        unset={unset}
        onBlur={undefined}
        onKeyDown={keepKeysInField}
        invalid={false}
        className={undefined}
      />
    </div>
  );
}

function renderNoEntries() {
  return <span className="text-muted-foreground px-1 text-xs">No entries</span>;
}

function keepKeysInField(event: React.KeyboardEvent) {
  if (event.key !== "Tab") {
    event.stopPropagation();
  }
}

interface RecordRow {
  id: number;
  record: SparseConfig;
}

function scalarText(value: JsonValue): string {
  return Array.isArray(value) || isJsonObject(value) ? JSON.stringify(value) : `${value}`;
}

function useSetting(spec: SettingSpec) {
  const { field, fieldState } = useController<FormConfig>({ name: spec.key });
  const watched = useWatch<FormConfig>({ name: spec.key });
  const value = watched;
  const { onChange } = field;
  const set = useCallback((next: JsonValue) => onChange(toFormValue(next)), [onChange]);
  const unset = useCallback(() => onChange(undefined), [onChange]);

  return { value, set, unset, onBlur: field.onBlur, invalid: fieldState.error !== undefined };
}
