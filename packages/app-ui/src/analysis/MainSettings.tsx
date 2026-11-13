import type {
  ActiveChoice,
  AppCommand,
  ChoiceName,
  ChoiceOptionName,
  EffectiveValue,
  JsonValue,
  MainSettingRow,
  RowControl,
  RowSwitch,
  RunCheck,
  SettingChoice,
  SettingPatch,
  SettingSpec,
  SparseConfig,
} from "@neherlab/app-contracts";
import { Badge } from "@neherlab/design/ui/badge";
import { Input } from "@neherlab/design/ui/input";
import { Label } from "@neherlab/design/ui/label";
import { Slider } from "@neherlab/design/ui/slider";
import { Switch } from "@neherlab/design/ui/switch";
import { cn } from "cn";
import { useCallback, useMemo, useState } from "react";
import { useFormContext } from "react-hook-form";

import { OptionToggle, type ToggleOption } from "../components/OptionToggle";
import { commandSettings } from "../settings/catalog";
import { choiceRow } from "../settings/choices";
import { documentedDefault, isChanged, resetValue, settingValue } from "../settings/config";
import { SettingHelp } from "../settings/help";
import { isNumber, isString, sameJson } from "../settings/json";
import { formatList, parseList } from "../settings/lists";
import { useDraftStore } from "../store/draft";
import { CheckItem } from "./ChecksPanel";
import { toFormValue, type FormConfig } from "./formValues";
import { NumberInput } from "./NumberInput";
import { SettingControl } from "./SettingControl";

const CLOCK_FILTER_MAX = 6;

const CLOCK_FILTER_STEP = 0.5;

const NO_CHECKS: readonly RunCheck[] = [];

const NO_EFFECTIVE: readonly EffectiveValue[] = [];

const EFFECTIVE_SHOWN: ReadonlySet<RowControl["kind"]> = new Set(["choice", "setting"]);

interface RowContext {
  command: AppCommand;
  config: SparseConfig;
  checks: readonly RunCheck[];
  effective: ReadonlyMap<string, string>;
  specs: ReadonlyMap<string, SettingSpec>;
  choices: readonly SettingChoice[];
  reported: readonly ActiveChoice[] | undefined;
  checked: SparseConfig | undefined;
  patch: (patches: readonly SettingPatch[]) => void;
  set: (key: string, value: JsonValue | undefined) => void;
  reset: (keys: readonly string[]) => void;
  get: (key: string) => JsonValue | undefined;
}

export function MainSettings({
  command,
  config,
  checks,
  reported,
  checked,
  effective,
}: {
  command: AppCommand;
  config: SparseConfig;
  checks: readonly RunCheck[] | undefined;
  reported: readonly ActiveChoice[] | undefined;
  checked: SparseConfig | undefined;
  effective: readonly EffectiveValue[] | undefined;
}) {
  const { setValue } = useFormContext<FormConfig>();
  const settings = commandSettings(command);

  const context = useMemo((): RowContext => {
    const specs = new Map(settings.settings.map((spec) => [spec.key, spec]));

    const write = (key: string, value: JsonValue | undefined) =>
      setValue(key, toFormValue(value), { shouldDirty: true, shouldValidate: true });

    const reset = (keys: readonly string[]) => {
      for (const key of keys) {
        const spec = specs.get(key);

        if (spec !== undefined) {
          write(key, resetValue(spec));
        }
      }
    };

    return {
      command,
      config,
      checks: checks ?? NO_CHECKS,
      effective: new Map((effective ?? NO_EFFECTIVE).map(({ key, text }) => [key, text])),
      specs,
      choices: settings.choices,
      reported,
      checked,
      patch: (patches) => {
        for (const { key, value } of patches) {
          write(key, value);
        }
      },
      reset,
      set: (key, value) => {
        const spec = specs.get(key);

        write(key, value);

        if (spec !== undefined && !sameJson(value, spec.default_value)) {
          reset(spec.conflicts);
        }
      },
      get: (key) => {
        const spec = specs.get(key);

        return spec === undefined ? undefined : settingValue(config, spec);
      },
    };
  }, [checked, checks, command, config, effective, reported, setValue, settings]);

  return (
    <div>
      {settings.main_rows.map((row) => (
        <Row key={row.label} context={context} row={row} />
      ))}
    </div>
  );
}

function Row({ context, row }: { context: RowContext; row: MainSettingRow }) {
  const { control } = row;

  switch (control.kind) {
    case "choice":
      return <ChoiceRow context={context} row={row} choice={control.choice} />;
    case "switches":
      return <SwitchesRow context={context} row={row} switches={control.switches} />;
    case "coalescent":
      return <CoalescentRow context={context} row={row} control={control} />;
    case "clock-filter":
      return <ClockFilterRow context={context} row={row} control={control} />;
    case "substitution-model":
      return <ModelRow context={context} row={row} control={control} />;
    case "setting":
      break;
  }

  return <SimpleRow context={context} row={row} />;
}

function SwitchesRow({
  context,
  row,
  switches,
}: {
  context: RowContext;
  row: MainSettingRow;
  switches: readonly RowSwitch[];
}) {
  return (
    <MainRow context={context} row={row}>
      <div className="grid gap-2">
        {switches.map(({ key, text, inverted }) =>
          inverted ? (
            <OnUnlessSwitch key={key} context={context} settingKey={key} text={text} />
          ) : (
            <SwitchLabel key={key} context={context} settingKey={key} text={text} />
          ),
        )}
      </div>
    </MainRow>
  );
}

function ChoiceRow({
  context,
  row: mainRow,
  choice: name,
}: {
  context: RowContext;
  row: MainSettingRow;
  choice: ChoiceName;
}) {
  const { choices, reported, checked, config, patch } = context;
  const choice = choices.find((candidate) => candidate.choice === name);
  const picked = useDraftStore((state) => state.draft.choices?.find((active) => active.choice === name)?.option);
  const pickChoice = useDraftStore((state) => state.pickChoice);

  const row = useMemo(
    () => (choice === undefined ? undefined : choiceRow(choice, { reported, picked, checked, config })),
    [checked, choice, config, picked, reported],
  );

  const options = useMemo(
    () =>
      (row?.options ?? []).map((option) => ({
        value: option,
        label: choice?.options.find((candidate) => candidate.option === option)?.label ?? option,
      })),
    [choice, row],
  );

  const onPick = useCallback(
    (option: ChoiceOptionName) => {
      pickChoice({ choice: name, option });
      patch(choice?.options.find((candidate) => candidate.option === option)?.patch ?? []);
    },
    [choice, name, patch, pickChoice],
  );

  if (choice === undefined || row?.selected === undefined) {
    return null;
  }

  const note = choice.options.find((candidate) => candidate.option === row.selected)?.note;

  return (
    <MainRow context={context} row={mainRow}>
      <OptionToggle label={mainRow.label} value={row.selected} onChange={onPick} options={options} />
      {row.shown.length > 0 && (
        <Inline>
          {row.shown.map((key) => (
            <Labeled key={key} text={rowLabel(context, mainRow, key)}>
              <Control context={context} settingKey={key} />
            </Labeled>
          ))}
        </Inline>
      )}
      {note !== undefined && <Note>{note}</Note>}
    </MainRow>
  );
}

function CoalescentRow({
  context,
  row,
  control,
}: {
  context: RowContext;
  row: MainSettingRow;
  control: Extract<RowControl, { kind: "coalescent" }>;
}) {
  const { get, set, reset } = context;
  const { setting, skyline_settings: skylineSettings, years_suffix: yearsSuffix, options } = control;
  const value = get(setting);
  const named = useMemo(() => namedCoalescents(context.specs.get(setting)), [context.specs, setting]);
  const derived = coalescentMode(value, named);
  const [picked, setPicked] = useState<string>();
  const mode = picked === "fixed" && derived === "none" ? "fixed" : derived;
  const years = yearsOf(value, yearsSuffix);

  const toggleOptions = useMemo(
    (): readonly ToggleOption<string>[] => options.map(({ value: option, label }) => ({ value: option, label })),
    [options],
  );

  const onPick = useCallback(
    (next: string) => {
      setPicked(next);

      if (next !== "skyline") {
        reset(skylineSettings);
      }

      if (named.has(next)) {
        set(setting, next);
      } else if (next === "none" || coalescentMode(get(setting), named) !== "fixed") {
        reset([setting]);
      }
    },
    [get, named, reset, set, setting, skylineSettings],
  );

  const onYears = useCallback(
    (next: number | string | undefined) => {
      if (next === undefined || next === "") {
        reset([setting]);
      } else {
        set(setting, `${next}${yearsSuffix}`);
      }
    },
    [reset, set, setting, yearsSuffix],
  );

  const note = options.find((option) => option.value === mode)?.note;

  return (
    <MainRow context={context} row={row}>
      <OptionToggle label={row.label} value={mode} onChange={onPick} options={toggleOptions} />
      {mode === "fixed" && (
        <>
          <Inline>
            <Labeled text={rowLabel(context, row, setting)}>
              <NumberInput
                aria-label={`${rowLabel(context, row, setting)} in years`}
                value={years}
                onValueChange={onYears}
                className="h-8 w-24 font-mono"
              />
              years
            </Labeled>
          </Inline>
          {years === undefined && value !== undefined && derived === "fixed" && (
            <Note>
              The configuration sets {isNumber(value) ? `${value} substitutions per site` : fixedText(value)}; a value
              in years replaces it.
            </Note>
          )}
        </>
      )}
      {mode === "skyline" && (
        <Inline>
          {skylineSettings.map((key) => (
            <Labeled key={key} text={rowLabel(context, row, key)}>
              <Control context={context} settingKey={key} className="w-20" />
            </Labeled>
          ))}
        </Inline>
      )}
      {note !== undefined && <Note>{note}</Note>}
    </MainRow>
  );
}

function namedCoalescents(spec: SettingSpec | undefined): ReadonlySet<string> {
  return new Set(spec?.options.map(({ value }) => value) ?? []);
}

function coalescentMode(value: JsonValue | undefined, named: ReadonlySet<string>): string {
  if (value === undefined || value === null) {
    return "none";
  }

  return isString(value) && named.has(value) ? value : "fixed";
}

function yearsOf(value: JsonValue | undefined, suffix: string): number | undefined {
  const number = isString(value) && value.endsWith(suffix) ? value.slice(0, -suffix.length) : "";
  const years = number.trim() === "" ? Number.NaN : Number(number);

  return Number.isFinite(years) && years > 0 ? years : undefined;
}

function ClockFilterRow({
  context,
  row,
  control,
}: {
  context: RowContext;
  row: MainSettingRow;
  control: Extract<RowControl, { kind: "clock-filter" }>;
}) {
  const { get, set } = context;
  const filter = get(control.setting);
  const value = isNumber(filter) ? filter : 0;

  const slid = useMemo(() => [Math.min(value, CLOCK_FILTER_MAX)], [value]);

  const onSlide = useCallback(
    (next: number | readonly number[]) => set(control.setting, [next].flat()[0] ?? 0),
    [control.setting, set],
  );

  return (
    <MainRow context={context} row={row}>
      <Inline>
        <Slider
          min={0}
          max={CLOCK_FILTER_MAX}
          step={CLOCK_FILTER_STEP}
          value={slid}
          aria-label={row.label}
          onValueChange={onSlide}
          className="w-56"
        />
        <Control context={context} settingKey={control.setting} className="w-20" />
      </Inline>
      <Note>{value === 0 ? control.off : control.on}</Note>
    </MainRow>
  );
}

function ModelRow({
  context,
  row,
  control,
}: {
  context: RowContext;
  row: MainSettingRow;
  control: Extract<RowControl, { kind: "substitution-model" }>;
}) {
  const { get, set, specs } = context;
  const { model_setting: modelKey, parameters_setting: paramsKey, parameters_example: example } = control;
  const modelSpec = specs.get(modelKey);
  const model = get(modelKey) ?? (modelSpec === undefined ? undefined : documentedDefault(modelSpec));
  const params = get(paramsKey);
  const modelHelp = modelSpec?.options.find((option) => option.value === model)?.help ?? "";
  const paramsLabel = rowLabel(context, row, paramsKey);

  const onParams = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => set(paramsKey, parseList(event.target.value, "string")),
    [paramsKey, set],
  );

  return (
    <MainRow context={context} row={row}>
      <Control context={context} settingKey={modelKey} className="w-auto" />
      {modelHelp !== "" && <Note>{modelHelp}</Note>}
      {specs.has(paramsKey) && (
        <Labeled text={paramsLabel}>
          <Input
            type="text"
            value={formatList(params)}
            aria-label={paramsLabel}
            placeholder={example}
            onChange={onParams}
            className="h-8 w-72 max-w-full font-mono"
          />
        </Labeled>
      )}
    </MainRow>
  );
}

function SimpleRow({ context, row }: { context: RowContext; row: MainSettingRow }) {
  const key = row.settings[0]?.key;
  const spec = key === undefined ? undefined : context.specs.get(key);

  if (key === undefined || spec === undefined) {
    return null;
  }

  return (
    <MainRow context={context} row={row}>
      <Control context={context} settingKey={key} className="max-w-72" />
      <SettingHelp spec={spec} />
    </MainRow>
  );
}

function MainRow({ context, row, children }: { context: RowContext; row: MainSettingRow; children: React.ReactNode }) {
  const keys = row.settings.map(({ key }) => key);

  const specs = keys.flatMap((key) => {
    const spec = context.specs.get(key);

    return spec === undefined ? [] : [spec];
  });

  const changed = specs.some((spec) => isChanged(context.config, spec));
  const checks = context.checks.filter((check) => check.settings.some((setting) => keys.includes(setting)));

  const effective =
    changed || !EFFECTIVE_SHOWN.has(row.control.kind)
      ? undefined
      : keys.map((key) => context.effective.get(key)).find((text) => text !== undefined);

  return (
    <div className="grid gap-2.5 border-t py-3 first:border-t-0 @2xl:grid-cols-[12.5rem_minmax(0,1fr)]">
      <div className="grid content-start gap-0.5">
        <span className="flex items-center gap-1.5 font-bold">
          {row.label}
          {changed && (
            <Badge variant="secondary" className="text-primary h-4 px-1.5 text-[0.6875rem]">
              changed
            </Badge>
          )}
        </span>
        {row.description !== undefined && <span className="text-muted-foreground text-xs">{row.description}</span>}
        {effective !== undefined && <span className="text-xs">{effective}</span>}
      </div>
      <div className="grid justify-items-start gap-2">
        {children}
        {checks.length > 0 && (
          <ul className="grid gap-1.5 text-xs">
            {checks.map((check) => (
              <CheckItem key={check.id} check={check} />
            ))}
          </ul>
        )}
      </div>
    </div>
  );
}

function Control({ context, settingKey, className }: { context: RowContext; settingKey: string; className?: string }) {
  const spec = context.specs.get(settingKey);

  return spec === undefined ? null : (
    <SettingControl
      command={context.command}
      spec={spec}
      label={spec.label}
      scope="main"
      className={cn("w-32", className)}
    />
  );
}

function OnUnlessSwitch({ context, settingKey, text }: { context: RowContext; settingKey: string; text: string }) {
  const { get, set, reset } = context;

  const onChange = useCallback(
    (checked: boolean) => (checked ? reset([settingKey]) : set(settingKey, true)),
    [reset, set, settingKey],
  );

  return (
    <Label className="font-normal">
      <Switch checked={get(settingKey) !== true} onCheckedChange={onChange} />
      {text}
    </Label>
  );
}

function SwitchLabel({ context, settingKey, text }: { context: RowContext; settingKey: string; text: string }) {
  const { get, set } = context;
  const onChange = useCallback((checked: boolean) => set(settingKey, checked), [set, settingKey]);

  return (
    <Label className="font-normal">
      <Switch checked={get(settingKey) === true} onCheckedChange={onChange} />
      {text}
    </Label>
  );
}

function Inline({ children }: { children: React.ReactNode }) {
  return <div className="flex flex-wrap items-center gap-2.5">{children}</div>;
}

function Labeled({ text, children }: { text: string; children: React.ReactNode }) {
  return (
    <Label className="text-muted-foreground max-w-full flex-wrap font-normal">
      {text}
      {children}
    </Label>
  );
}

function Note({ children }: { children: React.ReactNode }) {
  return <p className="text-muted-foreground text-xs">{children}</p>;
}

function fixedText(value: JsonValue): string {
  return isString(value) ? value : JSON.stringify(value);
}

function rowLabel(context: RowContext, row: MainSettingRow, key: string): string {
  return row.settings.find((setting) => setting.key === key)?.label ?? context.specs.get(key)?.label ?? key;
}
