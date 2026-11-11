import type {
  ActiveChoice,
  AppCommand,
  ChoiceName,
  ChoiceOptionName,
  JsonValue,
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
import type { SettingKey } from "../settings/commands";
import { documentedDefault, isChanged, resetValue, settingValue } from "../settings/config";
import { SettingHelp } from "../settings/help";
import { isNumber, isString, sameJson } from "../settings/json";
import { formatList, parseList } from "../settings/lists";
import { CheckItem } from "./ChecksPanel";
import { toFormValue, type FormConfig } from "./formValues";
import { NumberInput } from "./NumberInput";
import { SettingControl } from "./SettingControl";

type TimetreeKey = SettingKey<"timetree">;

const CLOCK_FILTER_MAX = 6;

const CLOCK_FILTER_STEP = 0.5;

const NO_CHECKS: readonly RunCheck[] = [];

const TOPOLOGY_SWITCHES = [
  { key: "no_merge_shared_mutations", text: "Merge polytomy siblings that share mutations" },
  { key: "no_hoist_reversions", text: "Hoist reversions in polytomies" },
  { key: "no_collapse_short_branches", text: "Collapse short internal branches" },
] as const satisfies readonly { key: SettingKey<"optimize">; text: string }[];

const TOPOLOGY_KEYS: ReadonlySet<string> = new Set(TOPOLOGY_SWITCHES.map(({ key }) => key));

const INTERVAL_KEYS = ["confidence"] as const satisfies readonly TimetreeKey[];

const CLOCK_FILTER_KEYS = ["clock_filter"] as const satisfies readonly TimetreeKey[];

const COALESCENT_KEYS = ["coalescent", "n_skyline", "skyline_stiffness"] as const satisfies readonly TimetreeKey[];

const POLYTOMY_KEYS = ["keep_polytomies"] as const satisfies readonly TimetreeKey[];

const MODEL_KEYS = ["model", "model_params"] as const satisfies readonly TimetreeKey[];

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
  root: {
    label: "Root",
    sub: "Pick rooting behavior",
    options: {
      reroot: { label: "Reroot", note: "Not set uses the command default." },
      outgroup: {
        label: "Outgroup",
        note: "The root is placed on the branch above the most recent common ancestor of the named tips.",
      },
      keep: {
        label: "Keep root",
        note: "The input root is kept; use this when the input tree is already rooted.",
      },
    },
    settings: { reroot: "Method", outgroup: "Tips" },
  },
};

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
  checks,
  reported,
  checked,
}: {
  command: AppCommand;
  config: SparseConfig;
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

        if (TOPOLOGY_KEYS.has(key)) {
          return settings.main_settings.find((first) => TOPOLOGY_KEYS.has(first)) === key ? (
            <TopologyRow key={key} context={context} shown={settings.main_settings} />
          ) : null;
        }

        return <SimpleRow key={key} context={context} settingKey={key} />;
      })}
    </div>
  );
}

function TimetreeSettings({ context }: { context: RowContext }) {
  return (
    <div className="py-1">
      <ChoiceRow context={context} choice="clock-rate" />
      <IntervalsRow context={context} />
      <CoalescentRow context={context} />
      <ChoiceRow context={context} choice="root" />
      <ClockFilterRow context={context} />
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
    <MainRow
      context={context}
      label="Confidence intervals"
      sub="Marginal dating uncertainty of nodes"
      keys={INTERVAL_KEYS}
    >
      <Inline>
        <SwitchLabel context={context} settingKey="confidence" text="Compute confidence intervals" />
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

function CoalescentRow({ context }: { context: RowContext }) {
  const { get, set, reset } = context;
  const value = get("coalescent");
  const derived = coalescentMode(value);
  const [picked, setPicked] = useState<string>();
  const mode = picked === "fixed" && derived === "none" ? "fixed" : derived;
  const years = yearsOf(value);

  const onPick = useCallback(
    (next: string) => {
      setPicked(next);

      if (next !== "skyline") {
        reset(SKYLINE_KEYS);
      }

      if (COALESCENT_NAMED.has(next)) {
        set("coalescent", next);
      } else if (next === "none" || coalescentMode(get("coalescent")) !== "fixed") {
        reset(["coalescent"]);
      }
    },
    [get, reset, set],
  );

  const onYears = useCallback(
    (next: number | string | undefined) => {
      if (next === undefined || next === "") {
        reset(["coalescent"]);
      } else {
        set("coalescent", `${next}y`);
      }
    },
    [reset, set],
  );

  return (
    <MainRow context={context} label="Coalescent model" sub="Tree prior on node times" keys={COALESCENT_KEYS}>
      <OptionToggle label="Coalescent model" value={mode} onChange={onPick} options={COALESCENT_OPTIONS} />
      {mode === "fixed" && (
        <>
          <Inline>
            <Labeled text="Time scale">
              <NumberInput
                aria-label="Coalescent time scale in years"
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
          <Labeled text="Skyline grid points">
            <Control context={context} settingKey="n_skyline" className="w-20" />
          </Labeled>
          <Labeled text="Skyline stiffness">
            <Control context={context} settingKey="skyline_stiffness" className="w-20" />
          </Labeled>
        </Inline>
      )}
    </MainRow>
  );
}

const COALESCENT_NAMED: ReadonlySet<string> = new Set(
  commandSettings("timetree")
    .settings.find(({ key }) => key === "coalescent")
    ?.options.map(({ value }) => value) ?? [],
);

const COALESCENT_LABELS = new Map([
  ["none", "None"],
  ["fixed", "Fixed constant"],
  ["opt", "Optimized constant"],
  ["skyline", "Skyline"],
]);

const COALESCENT_OPTIONS: readonly ToggleOption<string>[] = ["none", "fixed", ...COALESCENT_NAMED].map((value) => ({
  value,
  label: COALESCENT_LABELS.get(value) ?? value,
}));

const SKYLINE_KEYS = ["n_skyline", "skyline_stiffness"] as const satisfies readonly TimetreeKey[];

const YEARS_PATTERN = /^\+?(\d+(?:\.\d*)?|\.\d+)y$/u;

function coalescentMode(value: JsonValue | undefined): string {
  if (value === undefined || value === null) {
    return "none";
  }

  return isString(value) && COALESCENT_NAMED.has(value) ? value : "fixed";
}

function yearsOf(value: JsonValue | undefined): number | undefined {
  const match = isString(value) ? YEARS_PATTERN.exec(value) : null;

  return match?.[1] === undefined ? undefined : Number(match[1]);
}

function ClockFilterRow({ context }: { context: RowContext }) {
  const { get, set } = context;
  const filter = get("clock_filter");
  const value = isNumber(filter) ? filter : 0;

  const slid = useMemo(() => [Math.min(value, CLOCK_FILTER_MAX)], [value]);

  const onSlide = useCallback((next: number | readonly number[]) => set("clock_filter", [next].flat()[0] ?? 0), [set]);

  return (
    <MainRow
      context={context}
      label="Clock filter"
      sub="Outlier threshold in interquartile ranges"
      keys={CLOCK_FILTER_KEYS}
    >
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
          : `Samples more than ${value} interquartile ranges from the regression are ignored`}
      </Note>
    </MainRow>
  );
}

function ModelRow({ context }: { context: RowContext }) {
  const { get, set, specs } = context;
  const modelSpec = specs.get("model");
  const model = get("model") ?? (modelSpec === undefined ? undefined : documentedDefault(modelSpec));
  const params = get("model_params");
  const modelHelp = modelSpec?.options.find((option) => option.value === model)?.help ?? "";

  const onParams = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => set("model_params", parseList(event.target.value, "string")),
    [set],
  );

  return (
    <MainRow context={context} label="Substitution model" keys={MODEL_KEYS}>
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

function TopologyRow({ context, shown }: { context: RowContext; shown: readonly string[] }) {
  const switches = useMemo(() => TOPOLOGY_SWITCHES.filter(({ key }) => shown.includes(key)), [shown]);
  const keys = useMemo(() => switches.map(({ key }) => key), [switches]);

  return (
    <MainRow
      context={context}
      label="Topology optimizations"
      sub="On by default; switch off to keep the topology as it is"
      keys={keys}
    >
      <div className="grid gap-2">
        {switches.map(({ key, text }) => (
          <OnUnlessSwitch key={key} context={context} settingKey={key} text={text} />
        ))}
      </div>
    </MainRow>
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

interface ChoiceText {
  label: string;
  sub: string;
  options: Partial<Record<ChoiceOptionName, { label: string; note?: string }>>;
  settings: Partial<Record<string, string>>;
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
