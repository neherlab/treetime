import type {
  ActiveChoice,
  AppCommand,
  ChoiceName,
  ChoiceOptionName,
  InputFacts,
  JsonValue,
  RunCheck,
  SettingChoice,
  SettingPatch,
  SettingSpec,
  SparseConfig,
} from "@neherlab/app-contracts";
import { useCallback, useMemo, useState } from "react";
import { useFormContext } from "react-hook-form";

import { OptionToggle } from "../components/OptionToggle";
import { commandSettings } from "../settings/catalog";
import { choiceRow } from "../settings/choices";
import type { SettingKey } from "../settings/commands";
import { isChanged, resetValue, settingValue } from "../settings/config";
import { isNumber, isString, sameJson } from "../settings/json";
import { formatList, parseList } from "../settings/lists";
import { Badge } from "../ui/badge";
import { cn } from "../ui/cn";
import { Input } from "../ui/input";
import { Label } from "../ui/label";
import { NativeSelect, NativeSelectOption } from "../ui/native-select";
import { Slider } from "../ui/slider";
import { Switch } from "../ui/switch";
import { CheckItem } from "./ChecksPanel";
import { toFormValue, type FormConfig } from "./formValues";
import { NumberInput } from "./NumberInput";
import { SettingControl } from "./SettingControl";
import { SettingHelp } from "./SettingField";

type TimetreeKey = SettingKey<"timetree">;

const CLOCK_FILTER_MAX = 6;

const CLOCK_FILTER_STEP = 0.5;

const NO_CHECKS: readonly RunCheck[] = [];

const INTERVAL_KEYS = ["confidence", "covariation"] as const satisfies readonly TimetreeKey[];

const CLOCK_FILTER_KEYS = ["clock_filter"] as const satisfies readonly TimetreeKey[];

const RELAX_KEYS = ["relax"] as const satisfies readonly TimetreeKey[];

const POLYTOMY_KEYS = ["keep_polytomies"] as const satisfies readonly TimetreeKey[];

const MODEL_KEYS = ["model", "model_params"] as const satisfies readonly TimetreeKey[];

const ATTRIBUTE_KEYS = ["attribute"] as const satisfies readonly SettingKey<"mugration">[];

const CHOICE_TEXT: Record<ChoiceName, ChoiceText> = {
  "clock-rate": {
    label: "Clock rate",
    sub: "Substitutions per site per year",
    options: {
      estimate: {
        label: "Estimate from data",
        note: "Estimated by root-to-tip regression, then refined during dating.",
      },
      fixed: { label: "Fixed rate" },
    },
    settings: { clock_rate: "Rate", clock_std_dev: "Std. dev." },
  },
  "coalescent-prior": {
    label: "Coalescent prior",
    sub: "Prior on node times",
    options: {
      none: { label: "None", note: "Node times are set by the sequences and sampling dates only." },
      fixed: { label: "Fixed Tc" },
      optimized: { label: "Optimize Tc", note: "One constant Tc, fitted to the tree." },
      skyline: { label: "Skyline" },
    },
    settings: { coalescent: "Tc in years", skyline_n_points: "Grid points", skyline_stiffness: "Stiffness" },
  },
  root: {
    label: "Root",
    sub: "Where the tree is rooted",
    options: {
      reroot: { label: "Optimize root", note: "Not set uses the command default." },
      keep: {
        label: "Keep input root",
        note: "The input root is kept; use this when the tree is rooted with an outgroup.",
      },
    },
    settings: { reroot: "Method" },
  },
};

const NO_COLUMNS: readonly string[] = [];

interface RowContext {
  command: AppCommand;
  config: SparseConfig;
  checks: readonly RunCheck[];
  specs: ReadonlyMap<string, SettingSpec>;
  choices: readonly SettingChoice[];
  reported: readonly ActiveChoice[] | undefined;
  checked: SparseConfig | undefined;
  patch: (patches: readonly SettingPatch[]) => void;
  set: (key: string, value: JsonValue | undefined) => void;
  reset: (keys: readonly string[]) => void;
  get: (key: string) => JsonValue | undefined;
  example: (key: string) => JsonValue | undefined;
}

export function MainSettings({
  command,
  config,
  facts,
  checks,
  reported,
  checked,
}: {
  command: AppCommand;
  config: SparseConfig;
  facts: InputFacts | undefined;
  checks: readonly RunCheck[] | undefined;
  reported: readonly ActiveChoice[] | undefined;
  checked: SparseConfig | undefined;
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
      specs,
      choices: settings.choices,
      reported,
      checked,
      patch: (patches) => {
        for (const { path, value } of patches) {
          write(path.join("."), value);
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
      example: (key) => specs.get(key)?.examples[0],
    };
  }, [checked, checks, command, config, reported, setValue, settings]);

  if (command === "timetree") {
    return <TimetreeSettings context={context} />;
  }

  return (
    <div className="py-1">
      {settings.main_settings.map((key) => {
        const choice = settings.choices.find((candidate) => candidate.keys.includes(key));

        if (choice !== undefined) {
          return settings.main_settings.find((first) => choice.keys.includes(first)) === key ? (
            <ChoiceRow key={key} context={context} choice={choice.choice} />
          ) : null;
        }

        return key === "attribute" ? (
          <AttributeRow key={key} context={context} columns={facts?.metadata?.columns ?? NO_COLUMNS} />
        ) : (
          <SimpleRow key={key} context={context} settingKey={key} />
        );
      })}
    </div>
  );
}

function TimetreeSettings({ context }: { context: RowContext }) {
  return (
    <div className="py-1">
      <ChoiceRow context={context} choice="clock-rate" />
      <IntervalsRow context={context} />
      <ChoiceRow context={context} choice="coalescent-prior" />
      <ChoiceRow context={context} choice="root" />
      <ClockFilterRow context={context} />
      <RelaxRow context={context} />
      <MainRow context={context} label="Polytomies" sub="Nodes with more than two children" keys={POLYTOMY_KEYS}>
        <SwitchLabel context={context} settingKey="keep_polytomies" text="Keep polytomies unresolved" />
      </MainRow>
      <ModelRow context={context} />
      <SimpleRow context={context} settingKey="max_iter" />
    </div>
  );
}

function IntervalsRow({ context }: { context: RowContext }) {
  return (
    <MainRow context={context} label="Date intervals" sub="Marginal intervals on every node" keys={INTERVAL_KEYS}>
      <Inline>
        <SwitchLabel context={context} settingKey="confidence" text="Compute date intervals" />
        <SwitchLabel context={context} settingKey="covariation" text="Covariation-aware regression" />
      </Inline>
    </MainRow>
  );
}

function ChoiceRow({ context, choice: name }: { context: RowContext; choice: ChoiceName }) {
  const { choices, reported, checked, config, patch } = context;
  const choice = choices.find((candidate) => candidate.choice === name);
  const text = CHOICE_TEXT[name];
  const [picked, setPicked] = useState<ChoiceOptionName>();

  const row = useMemo(
    () => (choice === undefined ? undefined : choiceRow(choice, { reported, picked, checked, config })),
    [checked, choice, config, picked, reported],
  );

  const options = useMemo(
    () => (row?.options ?? []).map((option) => ({ value: option, label: text.options[option]?.label ?? option })),
    [row, text],
  );

  const onPick = useCallback(
    (option: ChoiceOptionName) => {
      setPicked(option);
      patch(choice?.options.find((candidate) => candidate.option === option)?.patch ?? []);
    },
    [choice, patch],
  );

  if (choice === undefined || row?.selected === undefined) {
    return null;
  }

  const note = text.options[row.selected]?.note;

  return (
    <MainRow context={context} label={text.label} sub={text.sub} keys={choice.keys}>
      <OptionToggle label={text.label} value={row.selected} onChange={onPick} options={options} />
      {row.shown.length > 0 && (
        <Inline>
          {row.shown.map((key) => (
            <Labeled key={key} text={text.settings[key] ?? context.specs.get(key)?.label ?? key}>
              <Control context={context} settingKey={key} />
            </Labeled>
          ))}
        </Inline>
      )}
      {note !== undefined && <Note>{note}</Note>}
    </MainRow>
  );
}

function ClockFilterRow({ context }: { context: RowContext }) {
  const { get, set } = context;
  const filter = get("clock_filter");
  const value = isNumber(filter) ? filter : 0;

  const slid = useMemo(() => [Math.min(value, CLOCK_FILTER_MAX)], [value]);

  const onSlide = useCallback((next: number | readonly number[]) => set("clock_filter", [next].flat()[0] ?? 0), [set]);

  return (
    <MainRow context={context} label="Clock filter" sub="Outlier threshold in IQD" keys={CLOCK_FILTER_KEYS}>
      <Inline>
        <Slider
          min={0}
          max={CLOCK_FILTER_MAX}
          step={CLOCK_FILTER_STEP}
          value={slid}
          aria-label="Clock filter threshold"
          onValueChange={onSlide}
          className="w-56"
        />
        <Control context={context} settingKey="clock_filter" className="w-20" />
      </Inline>
      <Note>
        {value === 0
          ? "Off: all samples stay in the regression"
          : `Samples more than ${value} interquartile distances from the regression are ignored`}
      </Note>
    </MainRow>
  );
}

function RelaxRow({ context }: { context: RowContext }) {
  const { get, set, reset, example, specs } = context;
  const relax = get("relax");
  const relaxed = Array.isArray(relax) && relax.length > 0;
  const valueNames = specs.get("relax")?.value_names;

  const [slack = "Slack", coupling = "Coupling"] = useMemo(() => (valueNames ?? []).map(valueNameLabel), [valueNames]);

  const onRelax = useCallback(
    (checked: boolean) => (checked ? set("relax", example("relax")) : reset(RELAX_KEYS)),
    [example, reset, set],
  );

  const onPair = useCallback((next: JsonValue) => set("relax", next), [set]);

  return (
    <MainRow context={context} label="Relaxed clock" sub="Rate variation between branches" keys={RELAX_KEYS}>
      <Inline>
        <Label className="font-normal">
          <Switch checked={relaxed} onCheckedChange={onRelax} />
          Relax the clock
        </Label>
        {relaxed && (
          <>
            <Labeled text={slack}>
              <NumberPairInput value={relax} index={0} onChange={onPair} label={slack} />
            </Labeled>
            <Labeled text={coupling}>
              <NumberPairInput value={relax} index={1} onChange={onPair} label={coupling} />
            </Labeled>
          </>
        )}
      </Inline>
      <Note>Values near 1 are weak priors. Coupling 0 is an uncorrelated clock.</Note>
    </MainRow>
  );
}

function ModelRow({ context }: { context: RowContext }) {
  const { get, set, specs } = context;
  const model = get("model");
  const params = get("model_params");
  const modelHelp = specs.get("model")?.options.find((option) => option.value === model)?.help ?? "";

  const onParams = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => set("model_params", parseList(event.target.value, "string")),
    [set],
  );

  return (
    <MainRow context={context} label="Substitution model" sub="Nucleotide or amino-acid model" keys={MODEL_KEYS}>
      <Control context={context} settingKey="model" className="w-auto" />
      {modelHelp !== "" && <Note>{modelHelp}</Note>}
      <Labeled text="Parameters">
        <Input
          type="text"
          value={formatList(params)}
          aria-label="Model parameters"
          placeholder="kappa=0.2 pis=0.25,0.25,0.25,0.25"
          onChange={onParams}
          className="h-8 w-72 max-w-full font-mono"
        />
      </Labeled>
    </MainRow>
  );
}

function AttributeRow({ context, columns }: { context: RowContext; columns: readonly string[] }) {
  const { get, set } = context;
  const value = get("attribute");

  const onColumn = useCallback(
    (event: React.ChangeEvent<HTMLSelectElement>) =>
      set("attribute", event.target.value === "" ? undefined : event.target.value),
    [set],
  );

  return (
    <MainRow context={context} label="Trait" sub="Metadata column to reconstruct" keys={ATTRIBUTE_KEYS}>
      {columns.length > 0 ? (
        <Labeled text="Column">
          <NativeSelect size="sm" aria-label="Trait column" value={isString(value) ? value : ""} onChange={onColumn}>
            <NativeSelectOption value="">Choose a column</NativeSelectOption>
            {columns.map((column) => (
              <NativeSelectOption key={column} value={column}>
                {column}
              </NativeSelectOption>
            ))}
          </NativeSelect>
        </Labeled>
      ) : (
        <Control context={context} settingKey="attribute" />
      )}
    </MainRow>
  );
}

function SimpleRow({ context, settingKey }: { context: RowContext; settingKey: string }) {
  const spec = context.specs.get(settingKey);
  const keys = useMemo(() => [settingKey], [settingKey]);

  if (spec === undefined) {
    return null;
  }

  return (
    <MainRow context={context} label={spec.label} keys={keys}>
      <Control context={context} settingKey={settingKey} className="max-w-72" />
      <SettingHelp spec={spec} />
    </MainRow>
  );
}

function MainRow({
  context,
  label,
  sub,
  keys,
  children,
}: {
  context: RowContext;
  label: string;
  sub?: string;
  keys: readonly string[];
  children: React.ReactNode;
}) {
  const specs = keys.flatMap((key) => {
    const spec = context.specs.get(key);

    return spec === undefined ? [] : [spec];
  });

  const changed = specs.some((spec) => isChanged(context.config, spec));
  const checks = context.checks.filter((check) => check.settings.some((setting) => keys.includes(setting)));

  return (
    <div className="grid gap-3.5 border-t px-3.5 py-3 first:border-t-0 @xl:grid-cols-[12.5rem_minmax(0,1fr)]">
      <div className="grid content-start gap-0.5">
        <span className="flex items-center gap-1.5 font-bold">
          {label}
          {changed && (
            <Badge variant="secondary" className="text-primary h-4 px-1.5 text-[0.6875rem]">
              changed
            </Badge>
          )}
        </span>
        {sub !== undefined && <span className="text-muted-foreground text-xs">{sub}</span>}
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
        <div className="flex flex-wrap gap-1">
          {specs.map((spec) => (
            <Badge key={spec.key} variant="outline" className="text-muted-foreground h-4 px-1 font-mono font-normal">
              {spec.flag}
            </Badge>
          ))}
        </div>
      </div>
    </div>
  );
}

function Control({ context, settingKey, className }: { context: RowContext; settingKey: string; className?: string }) {
  const spec = context.specs.get(settingKey);

  return spec === undefined ? null : (
    <SettingControl command={context.command} spec={spec} label={spec.label} className={cn("w-32", className)} />
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

function NumberPairInput({
  value,
  index,
  onChange,
  label,
}: {
  value: JsonValue | undefined;
  index: number;
  onChange: (value: JsonValue) => void;
  label: string;
}) {
  const pair = useMemo(() => (Array.isArray(value) ? value : []), [value]);

  const onInput = useCallback(
    (next: number | string | undefined) =>
      onChange(pair.map((current, position) => (position === index ? (next ?? "") : current))),
    [index, onChange, pair],
  );

  return <NumberInput aria-label={label} value={pair[index]} onValueChange={onInput} className="h-8 w-24 font-mono" />;
}

interface ChoiceText {
  label: string;
  sub: string;
  options: Partial<Record<ChoiceOptionName, { label: string; note?: string }>>;
  settings: Partial<Record<string, string>>;
}

function valueNameLabel(name: string): string {
  const lower = name.toLowerCase().replaceAll("_", " ");

  return lower.charAt(0).toUpperCase() + lower.slice(1);
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
