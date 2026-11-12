import type { JsonValue, SettingOption, SettingSpec } from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "@neherlab/design/ui/collapsible";
import { Popover, PopoverContent, PopoverTrigger } from "@neherlab/design/ui/popover";
import type { ReactNode } from "react";
import Info from "~icons/lucide/info";

import { m, useLocale, type Locale } from "../i18n/locale";
import { type MarkdownComponents, MarkdownText } from "../markdown";
import { documentedDefault } from "./config";
import { isJsonObject } from "./json";
import { formatList } from "./lists";

const CODE_CLASS = "bg-muted text-foreground rounded-sm px-1 py-px font-mono text-[0.92em]";

const INLINE_COMPONENTS: MarkdownComponents = {
  p: ({ children }) => <span>{children}</span>,
  code: ({ children }) => <code className={CODE_CLASS}>{children}</code>,
};

const BLOCK_COMPONENTS: MarkdownComponents = {
  p: ({ children }) => <p className="mt-1.5">{children}</p>,
  code: ({ children }) => <code className={CODE_CLASS}>{children}</code>,
};

export function SettingSummary({ spec, label }: { spec: SettingSpec; label: ReactNode }) {
  const locale = useLocale();

  return (
    <div className="flex flex-wrap items-baseline gap-2">
      {label} <code className="text-muted-foreground font-mono text-xs">{spec.flag}</code>{" "}
      {spec.role === "setting" && (
        <span className="text-muted-foreground text-xs">
          {m.settings_default({ value: defaultText(documentedDefault(spec), locale, spec.options) }, { locale })}
        </span>
      )}
    </div>
  );
}

export function SettingHelpPopover({ spec, value }: { spec: SettingSpec; value: JsonValue | undefined }) {
  const locale = useLocale();

  return (
    <Popover>
      <PopoverTrigger
        render={
          <Button
            type="button"
            variant="ghost"
            size="icon-xs"
            aria-label={m.settings_about({ label: spec.label }, { locale })}
            className="text-muted-foreground"
          />
        }
      >
        <Info aria-hidden />
      </PopoverTrigger>
      <PopoverContent align="start" className="grid w-96 max-w-[calc(100vw-2rem)] gap-2.5">
        <span className="font-bold">{spec.label}</span>
        <div className="text-muted-foreground grid gap-1.5">
          <MarkdownText
            text={[spec.help, spec.more].filter((text) => text !== "").join("\n\n")}
            components={BLOCK_COMPONENTS}
          />
        </div>
        {spec.options.length > 0 && (
          <dl className="grid grid-cols-[auto_minmax(0,1fr)] gap-x-3 gap-y-1 text-xs">
            {spec.options.map((option) => (
              <div key={option.value} className="contents">
                <dt>
                  {option.label} <code className={CODE_CLASS}>{option.value}</code>
                </dt>
                <dd className="text-muted-foreground">{option.help}</dd>
              </div>
            ))}
          </dl>
        )}
        <p className="text-xs">
          {m.settings_command_line({}, { locale })}: <code className={CODE_CLASS}>{commandLineText(spec, value)}</code>
        </p>
      </PopoverContent>
    </Popover>
  );
}

export function commandLineText(spec: SettingSpec, value: JsonValue | undefined): string {
  if (value === undefined || value === null || isJsonObject(value) || spec.kind === "records") {
    return [spec.flag, ...spec.value_names].join(" ");
  }

  if (value === true || value === false) {
    return spec.flag;
  }

  const values = Array.isArray(value) ? value : [value];

  return values.length === 0
    ? spec.flag
    : [
        spec.flag,
        ...values.map((item) => (isJsonObject(item) || Array.isArray(item) ? JSON.stringify(item) : String(item))),
      ].join(" ");
}

export function SettingHelpText({ spec }: { spec: SettingSpec }) {
  return (
    <div className="text-muted-foreground max-w-[72ch] text-[0.8125rem]">
      <MarkdownText text={spec.help} components={INLINE_COMPONENTS} />
    </div>
  );
}

export function SettingHelp({ spec }: { spec: SettingSpec }) {
  const locale = useLocale();

  return (
    <div className="text-muted-foreground max-w-[72ch] text-[0.8125rem]">
      <MarkdownText text={spec.help} components={INLINE_COMPONENTS} />
      {spec.more !== "" && (
        <Collapsible className="inline">
          <CollapsibleTrigger className="text-primary ml-1 underline-offset-4 hover:underline">
            {m.settings_more({}, { locale })}
          </CollapsibleTrigger>
          <CollapsibleContent hiddenUntilFound>
            <MarkdownText text={spec.more} components={BLOCK_COMPONENTS} />
          </CollapsibleContent>
        </Collapsible>
      )}
    </div>
  );
}

export function defaultText(
  value: JsonValue | undefined,
  locale: Locale,
  options: readonly SettingOption[] = [],
): string {
  if (value === undefined) {
    return m.settings_not_set({}, { locale });
  }

  if (Array.isArray(value)) {
    return value.length === 0
      ? m.settings_none({}, { locale })
      : formatList(value.map((item) => optionLabel(item, options)));
  }

  if (value === true) {
    return m.settings_on({}, { locale });
  }

  if (value === false) {
    return m.settings_off({}, { locale });
  }

  return isJsonObject(value) ? JSON.stringify(value) : String(optionLabel(value, options));
}

function optionLabel<T extends JsonValue>(value: T, options: readonly SettingOption[]): T | string {
  return options.find((option) => option.value === value)?.label ?? value;
}
