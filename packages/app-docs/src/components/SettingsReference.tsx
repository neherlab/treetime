import type { CommandSettings, InputNeed, JsonValue, SettingSpec } from "@neherlab/app-contracts";
import { fuzzyFilter } from "@neherlab/app-ui/fuzzy";
import { LocaleProvider, m, useLocale, type Locale } from "@neherlab/app-ui/i18n";
import { groupedSpecs } from "@neherlab/app-ui/settings/groups";
import { SettingHelp, SettingSummary } from "@neherlab/app-ui/settings/help";
import { Input } from "@neherlab/app-ui/ui/input";
import { Label } from "@neherlab/app-ui/ui/label";
import { type ChangeEvent, useCallback, useId, useMemo, useState } from "react";

export function SettingsReference({ settings, locale }: { settings: CommandSettings; locale: Locale }) {
  const filterId = useId();
  const [query, setQuery] = useState("");
  const onChange = useCallback((event: ChangeEvent<HTMLInputElement>) => setQuery(event.currentTarget.value), []);
  const matching = useMemo(() => matchingSpecs(settings.settings, query), [query, settings.settings]);
  const groups = useMemo(() => groupedSpecs(settings, matching), [matching, settings]);
  const flags = useMemo(() => new Map(settings.settings.map((spec) => [spec.key, spec.flag])), [settings.settings]);

  return (
    <LocaleProvider locale={locale}>
      <div className="not-prose grid gap-8">
        <CommandOverview settings={settings} />
        <div className="grid gap-2">
          <Label htmlFor={filterId}>{m.settings_filter_label({}, { locale })}</Label>
          <Input
            id={filterId}
            type="search"
            value={query}
            onChange={onChange}
            placeholder={m.settings_filter_placeholder({}, { locale })}
            className="max-w-md"
          />
          <output aria-live="polite" className="text-muted-foreground text-sm">
            {m.settings_filter_count({ count: matching.length }, { locale })}
          </output>
        </div>
        {groups.map(([group, specs]) => (
          <section key={group} aria-labelledby={groupId(group)} className="grid gap-4">
            <h2 id={groupId(group)} className="font-heading border-b pb-1 text-xl font-semibold">
              {group}
            </h2>
            {specs.map((spec) => (
              <SettingEntry key={spec.key} spec={spec} flags={flags} />
            ))}
          </section>
        ))}
      </div>
    </LocaleProvider>
  );
}

function CommandOverview({ settings }: { settings: CommandSettings }) {
  const locale = useLocale();

  return (
    <div className="grid gap-4">
      <p className="text-lg">{settings.description}</p>
      {settings.inputs.length > 0 && (
        <div className="grid gap-2">
          <h2 className="font-heading text-xl font-semibold">{m.settings_inputs({}, { locale })}</h2>
          <ul className="grid gap-1 text-sm">
            {settings.inputs.map((input) => (
              <li key={input.kind}>
                <span className="font-bold">{input.label}</span>{" "}
                <span className="text-muted-foreground">
                  ({needText(input.need, locale)}): {input.formats}
                </span>
              </li>
            ))}
          </ul>
        </div>
      )}
      <nav aria-label={m.settings_groups({}, { locale })}>
        <ul className="flex flex-wrap gap-x-4 gap-y-1 text-sm">
          {settings.groups.map((group) => (
            <li key={group}>
              <a href={`#${groupId(group)}`} className="text-primary underline-offset-3 hover:underline">
                {group}
              </a>
            </li>
          ))}
        </ul>
      </nav>
    </div>
  );
}

function SettingEntry({ spec, flags }: { spec: SettingSpec; flags: ReadonlyMap<string, string> }) {
  const locale = useLocale();

  return (
    <article id={spec.key} className="bg-card grid scroll-mt-20 gap-2 rounded-lg border p-4">
      <SettingSummary
        spec={spec}
        label={
          <h3 className="font-heading text-base font-semibold">
            <a href={`#${spec.key}`} className="hover:underline">
              {spec.label}
            </a>
          </h3>
        }
      />
      <dl className="grid grid-cols-[max-content_1fr] gap-x-4 gap-y-1 text-sm">
        {spec.value_names.length > 0 && (
          <Term label={m.settings_value({}, { locale })}>
            <code className="font-mono text-xs">{spec.value_names.join(" ")}</code>
          </Term>
        )}
        <Term label={m.settings_config_key({}, { locale })}>
          <code className="font-mono text-xs">{spec.key}</code>
        </Term>
        {spec.options.length > 0 && (
          <Term label={m.settings_options({}, { locale })}>
            <ul className="grid gap-0.5">
              {spec.options.map((option) => (
                <li key={option.value}>
                  <code className="font-mono text-xs">{option.value}</code>
                  {option.help !== "" && <span className="text-muted-foreground"> {option.help}</span>}
                </li>
              ))}
            </ul>
          </Term>
        )}
        {spec.conflicts.length > 0 && (
          <Term label={m.settings_conflicts({}, { locale })}>
            {spec.conflicts.map((key) => (
              <a key={key} href={`#${key}`} className="mr-2 font-mono text-xs underline-offset-3 hover:underline">
                {flags.get(key) ?? key}
              </a>
            ))}
          </Term>
        )}
        {spec.examples.length > 0 && (
          <Term label={m.settings_examples({}, { locale })}>
            <ul className="grid gap-0.5">
              {spec.examples.map((example) => (
                <li key={exampleText(example)}>
                  <code className="font-mono text-xs break-all">{exampleText(example)}</code>
                </li>
              ))}
            </ul>
          </Term>
        )}
      </dl>
      <SettingHelp spec={spec} />
    </article>
  );
}

function Term({ label, children }: { label: string; children: React.ReactNode }) {
  return (
    <>
      <dt className="text-muted-foreground">{label}</dt>
      <dd className="min-w-0">{children}</dd>
    </>
  );
}

const NEED_MESSAGES: Record<InputNeed, (locale: Locale) => string> = {
  required: (locale: Locale) => m.settings_need_required({}, { locale }),
  recommended: (locale: Locale) => m.settings_need_recommended({}, { locale }),
  optional: (locale: Locale) => m.settings_need_optional({}, { locale }),
};

function needText(need: InputNeed, locale: Locale): string {
  return NEED_MESSAGES[need](locale);
}

function matchingSpecs(specs: readonly SettingSpec[], query: string): SettingSpec[] {
  return query.trim() === "" ? [...specs] : fuzzyFilter(specs, query, specText);
}

function specText(spec: SettingSpec): string {
  return [spec.label, spec.flag, spec.key, spec.help].join(" ");
}

function exampleText(example: JsonValue): string {
  return JSON.stringify(example);
}

function groupId(group: string): string {
  return `group-${group.toLowerCase().replaceAll(/[^a-z0-9]+/gu, "-")}`;
}
