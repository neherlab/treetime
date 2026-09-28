import type { AppCommand, InputFacts, RunCheck } from "@neherlab/app-contracts";
import { useCallback, useMemo } from "react";
import { useFormContext } from "react-hook-form";

import { COMMAND_SETTINGS, type SettingSpec } from "../settings/catalog";
import { MAIN_SETTING_KEYS, type SettingKey } from "../settings/commands";
import { isChanged, resetValue, settingValue } from "../settings/config";
import { isNumber, isString, sameJson, type JsonObject, type JsonValue } from "../settings/json";
import { formatList, parseList } from "../settings/lists";
import { Segmented, Switch, cn, type SegmentedOption } from "../ui";
import { CheckItem } from "./ChecksPanel";
import { toFormValue, type FormConfig } from "./formValues";
import { SettingControl } from "./SettingControl";
import { SettingHelp } from "./SettingField";

type TimetreeKey = SettingKey<"timetree">;

const CLOCK_FILTER_MAX = 6;

const CLOCK_FILTER_STEP = 0.5;

const NO_CHECKS: readonly RunCheck[] = [];

type RateMode = "estimate" | "fixed";

type CoalescentMode = "none" | "fixed" | "optimized" | "skyline";

type RootMode = "optimize" | "keep";

const RATE_MODES: ReadonlyArray<SegmentedOption<RateMode>> = [
  { value: "estimate", label: "Estimate from data" },
  { value: "fixed", label: "Fixed rate" },
];

const COALESCENT_MODES: ReadonlyArray<SegmentedOption<CoalescentMode>> = [
  { value: "none", label: "None" },
  { value: "fixed", label: "Fixed Tc" },
  { value: "optimized", label: "Optimize Tc" },
  { value: "skyline", label: "Skyline" },
];

const ROOT_MODES: ReadonlyArray<SegmentedOption<RootMode>> = [
  { value: "optimize", label: "Optimize root" },
  { value: "keep", label: "Keep input root" },
];

const CLOCK_RATE_KEYS = ["clock_rate", "clock_std_dev"] as const satisfies readonly TimetreeKey[];

const INTERVAL_KEYS = ["confidence", "covariation"] as const satisfies readonly TimetreeKey[];

const COALESCENT_KEYS = [
  "coalescent",
  "coalescent_opt",
  "coalescent_skyline",
  "skyline_n_points",
  "skyline_stiffness",
] as const satisfies readonly TimetreeKey[];

const ROOT_KEYS = ["reroot", "keep_root"] as const satisfies readonly TimetreeKey[];

const CLOCK_FILTER_KEYS = ["clock_filter"] as const satisfies readonly TimetreeKey[];

const RELAX_KEYS = ["relax"] as const satisfies readonly TimetreeKey[];

const POLYTOMY_KEYS = ["keep_polytomies"] as const satisfies readonly TimetreeKey[];

const MODEL_KEYS = ["model", "model_params"] as const satisfies readonly TimetreeKey[];

const ATTRIBUTE_KEYS = ["attribute"] as const satisfies readonly SettingKey<"mugration">[];

const NO_COLUMNS: readonly string[] = [];

interface RowContext {
  command: AppCommand;
  config: JsonObject;
  checks: readonly RunCheck[];
  specs: ReadonlyMap<string, SettingSpec>;
  set: (key: string, value: JsonValue) => void;
  reset: (keys: readonly string[]) => void;
  get: (key: string) => JsonValue;
  example: (key: string) => JsonValue;
}

export function MainSettings({
  command,
  config,
  facts,
  checks,
}: {
  command: AppCommand;
  config: JsonObject;
  facts: InputFacts | undefined;
  checks: readonly RunCheck[] | undefined;
}) {
  const { setValue } = useFormContext<FormConfig>();

  const context = useMemo((): RowContext => {
    const specs = new Map(COMMAND_SETTINGS[command].specs.map((spec) => [spec.key, spec]));

    const write = (key: string, value: JsonValue) =>
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

        return spec === undefined ? null : settingValue(config, spec);
      },
      example: (key) => specs.get(key)?.examples[0] ?? null,
    };
  }, [checks, command, config, setValue]);

  if (command === "timetree") {
    return <TimetreeSettings context={context} />;
  }

  const keys: readonly string[] = MAIN_SETTING_KEYS[command];
  const rootKeys: readonly string[] = ROOT_KEYS;
  const hasRoot = keys.includes("reroot") && context.specs.has("keep_root");
  const rest = keys.filter((key) => !(hasRoot && rootKeys.includes(key)));

  return (
    <div className="py-1">
      {hasRoot && <RootRow context={context} />}
      {rest.map((key) =>
        key === "attribute" ? (
          <AttributeRow key={key} context={context} columns={facts?.metadata?.columns ?? NO_COLUMNS} />
        ) : (
          <SimpleRow key={key} context={context} settingKey={key} />
        ),
      )}
    </div>
  );
}

function TimetreeSettings({ context }: { context: RowContext }) {
  return (
    <div className="py-1">
      <ClockRateRow context={context} />
      <IntervalsRow context={context} />
      <CoalescentRow context={context} />
      <RootRow context={context} />
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

function ClockRateRow({ context }: { context: RowContext }) {
  const { get, set, reset, example } = context;
  const fixed = get("clock_rate") !== null;

  const onMode = useCallback(
    (mode: RateMode) => {
      reset(CLOCK_RATE_KEYS);

      if (mode === "fixed") {
        set("clock_rate", example("clock_rate"));
      }
    },
    [example, reset, set],
  );

  return (
    <MainRow context={context} label="Clock rate" sub="Substitutions per site per year" keys={CLOCK_RATE_KEYS}>
      <Segmented label="Clock rate" value={fixed ? "fixed" : "estimate"} onChange={onMode} options={RATE_MODES} />
      {fixed ? (
        <Inline>
          <Labeled text="Rate">
            <Control context={context} settingKey="clock_rate" />
          </Labeled>
          <Labeled text="Std. dev.">
            <Control context={context} settingKey="clock_std_dev" />
          </Labeled>
        </Inline>
      ) : (
        <Note>Estimated by root-to-tip regression, then refined during dating.</Note>
      )}
    </MainRow>
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

function CoalescentRow({ context }: { context: RowContext }) {
  const { set, reset, example } = context;
  const mode = coalescentMode(context);

  const onMode = useCallback(
    (next: CoalescentMode) => {
      reset(COALESCENT_KEYS);

      if (next === "fixed") {
        set("coalescent", example("coalescent"));
      } else if (next === "optimized") {
        set("coalescent_opt", true);
      } else if (next === "skyline") {
        set("coalescent_skyline", true);
      }
    },
    [example, reset, set],
  );

  return (
    <MainRow context={context} label="Coalescent prior" sub="Prior on node times" keys={COALESCENT_KEYS}>
      <Segmented label="Coalescent prior" value={mode} onChange={onMode} options={COALESCENT_MODES} />
      {mode === "fixed" && (
        <Inline>
          <Labeled text="Tc in years">
            <Control context={context} settingKey="coalescent" />
          </Labeled>
        </Inline>
      )}
      {mode === "skyline" && (
        <Inline>
          <Labeled text="Grid points">
            <Control context={context} settingKey="skyline_n_points" />
          </Labeled>
          <Labeled text="Stiffness">
            <Control context={context} settingKey="skyline_stiffness" />
          </Labeled>
        </Inline>
      )}
      {mode === "optimized" && <Note>One constant Tc, fitted to the tree.</Note>}
      {mode === "none" && <Note>Node times are set by the sequences and sampling dates only.</Note>}
    </MainRow>
  );
}

function coalescentMode(context: RowContext): CoalescentMode {
  if (context.get("coalescent_skyline") === true) {
    return "skyline";
  }

  if (context.get("coalescent_opt") === true) {
    return "optimized";
  }

  return context.get("coalescent") === null ? "none" : "fixed";
}

function RootRow({ context }: { context: RowContext }) {
  const { get, set } = context;
  const keep = get("keep_root") === true;

  const onMode = useCallback((mode: RootMode) => set("keep_root", mode === "keep"), [set]);

  return (
    <MainRow context={context} label="Root" sub="Where the tree is rooted" keys={ROOT_KEYS}>
      <Segmented label="Root" value={keep ? "keep" : "optimize"} onChange={onMode} options={ROOT_MODES} />
      {keep ? (
        <Note>The input root is kept; use this when the tree is rooted with an outgroup.</Note>
      ) : (
        <Inline>
          <Labeled text="Method">
            <Control context={context} settingKey="reroot" />
          </Labeled>
          <Note>Not set uses the command default.</Note>
        </Inline>
      )}
    </MainRow>
  );
}

function ClockFilterRow({ context }: { context: RowContext }) {
  const { get, set } = context;
  const value = Number(get("clock_filter") ?? 0);

  const onSlide = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => set("clock_filter", Number(event.target.value)),
    [set],
  );

  return (
    <MainRow context={context} label="Clock filter" sub="Outlier threshold in IQD" keys={CLOCK_FILTER_KEYS}>
      <Inline>
        <input
          type="range"
          min={0}
          max={CLOCK_FILTER_MAX}
          step={CLOCK_FILTER_STEP}
          value={Math.min(value, CLOCK_FILTER_MAX)}
          aria-label="Clock filter threshold"
          onChange={onSlide}
          className="accent-accent w-56"
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
        <span className="text-ink-muted inline-flex items-center gap-1.5">
          <Switch label="Relax the clock" checked={relaxed} onCheckedChange={onRelax} />
          Relax the clock
        </span>
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
        <input
          type="text"
          value={formatList(params)}
          aria-label="Model parameters"
          placeholder="kappa=0.2 pis=0.25,0.25,0.25,0.25"
          onChange={onParams}
          className="border-line-strong bg-surface-1 w-72 max-w-full rounded-md border px-2 py-1"
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
      set("attribute", event.target.value === "" ? null : event.target.value),
    [set],
  );

  return (
    <MainRow context={context} label="Trait" sub="Metadata column to reconstruct" keys={ATTRIBUTE_KEYS}>
      {columns.length > 0 ? (
        <Labeled text="Column">
          <select
            aria-label="Trait column"
            value={isString(value) ? value : ""}
            onChange={onColumn}
            className="border-line-strong bg-surface-1 rounded-md border px-2 py-1"
          >
            <option value="">Choose a column</option>
            {columns.map((column) => (
              <option key={column} value={column}>
                {column}
              </option>
            ))}
          </select>
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
    <div className="border-surface-3 grid gap-3.5 border-t px-3.5 py-3 first:border-t-0 @xl:grid-cols-[12.5rem_minmax(0,1fr)]">
      <div className="font-bold">
        {label}
        {changed && (
          <span
            title="Changed from default"
            className="bg-accent ml-1.5 inline-block size-1.5 rounded-full align-middle"
          />
        )}
        {sub !== undefined && <small className="text-ink-faint block text-xs font-normal">{sub}</small>}
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
            <code key={spec.key} className="bg-surface-3 text-ink-muted rounded-sm px-1 font-mono text-[0.6875rem]">
              {spec.flag}
            </code>
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
    <span className="text-ink-muted inline-flex items-center gap-1.5">
      <Switch label={text} checked={get(settingKey) === true} onCheckedChange={onChange} />
      {text}
    </span>
  );
}

function NumberPairInput({
  value,
  index,
  onChange,
  label,
}: {
  value: JsonValue;
  index: number;
  onChange: (value: JsonValue) => void;
  label: string;
}) {
  const pair = useMemo(() => (Array.isArray(value) ? value : []), [value]);
  const item = pair[index];

  const onInput = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) =>
      onChange(
        pair.map((current, position) => (position === index ? numberOrNull(event.target.valueAsNumber) : current)),
      ),
    [index, onChange, pair],
  );

  return (
    <input
      type="number"
      step="any"
      aria-label={label}
      value={isNumber(item) ? item : ""}
      onChange={onInput}
      className="border-line-strong bg-surface-1 w-24 rounded-md border px-2 py-1"
    />
  );
}

function numberOrNull(value: number): number | null {
  return Number.isNaN(value) ? null : value;
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
    <span className="text-ink-muted inline-flex max-w-full flex-wrap items-center gap-1.5">
      {text}
      {children}
    </span>
  );
}

function Note({ children }: { children: React.ReactNode }) {
  return <p className="text-ink-faint text-xs">{children}</p>;
}
