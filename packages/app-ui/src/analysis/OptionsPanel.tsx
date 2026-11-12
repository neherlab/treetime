import type {
  ActiveChoice,
  AppCommand,
  EffectiveValue,
  InputFacts,
  RunCheck,
  SparseConfig,
} from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { HoverCard, HoverCardContent, HoverCardTrigger } from "@neherlab/design/ui/hover-card";
import { Popover, PopoverContent, PopoverTrigger } from "@neherlab/design/ui/popover";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@neherlab/design/ui/select";
import { Sheet, SheetContent, SheetDescription, SheetHeader, SheetTitle } from "@neherlab/design/ui/sheet";
import { useCallback, useMemo, type ReactNode } from "react";
import Info from "~icons/lucide/info";
import SlidersHorizontal from "~icons/lucide/sliders-horizontal";

import { COMMANDS, commandSettings } from "../settings/catalog";
import { carryOverDraft } from "../settings/commands";
import { changedSpecs } from "../settings/config";
import { SettingHelp } from "../settings/help";
import { useDraftStore } from "../store/draft";
import { type DetectedFact, detectedFacts } from "./detectedFacts";
import { InputsSection } from "./InputsSection";
import { MainSettings } from "./MainSettings";
import { SettingControl } from "./SettingControl";
import { useSettingFocus } from "./settingFocus";
import { AllSettingsPanel } from "./SettingsPanel";

const ANALYSES = COMMANDS.map(({ command, title }) => ({ value: command, label: title }));

export function OptionsPanel({
  command,
  config,
  checks,
  reported,
  checked,
  facts,
  effective,
}: {
  command: AppCommand;
  config: SparseConfig;
  checks: readonly RunCheck[] | undefined;
  reported: readonly ActiveChoice[] | undefined;
  checked: SparseConfig | undefined;
  facts: InputFacts | undefined;
  effective: readonly EffectiveValue[] | undefined;
}) {
  const changed = changedSpecs(commandSettings(command).settings, config).length;

  useSettingFocus();

  return (
    <div className="grid gap-5 px-4 pt-3 pb-5">
      <AnalysisSelector command={command} config={config} />
      <Section title="Data">
        <InputsSection command={command} facts={facts} />
        {facts !== undefined && <DetectedFacts command={command} facts={detectedFacts(facts.detected)} />}
      </Section>
      <Section title="Main settings" hint={changed === 0 ? "All defaults" : `${changed} changed`}>
        <MainSettings
          command={command}
          config={config}
          checks={checks}
          reported={reported}
          checked={checked}
          effective={effective}
        />
        <AllSettings command={command} config={config} />
      </Section>
    </div>
  );
}

function AnalysisSelector({ command, config }: { command: AppCommand; config: SparseConfig }) {
  const settings = commandSettings(command);

  const select = useCallback(
    (value: string | null) => {
      const target = COMMANDS.find((candidate) => candidate.command === value)?.command;

      if (target === undefined || target === command) {
        return;
      }

      const { draft, load } = useDraftStore.getState();

      load({ ...draft, command: target, ...carryOverDraft(target, command, config, draft.sources) });
    },
    [command, config],
  );

  return (
    <div className="flex items-center gap-2">
      <Select items={ANALYSES} value={command} onValueChange={select}>
        <SelectTrigger aria-label="Analysis" className="min-w-0 flex-1 font-bold">
          <SelectValue />
        </SelectTrigger>
        <SelectContent>
          {ANALYSES.map(({ value, label }) => (
            <SelectItem key={value} value={value}>
              {label}
            </SelectItem>
          ))}
        </SelectContent>
      </Select>
      <HoverCard>
        <HoverCardTrigger
          render={<Button type="button" variant="ghost" size="icon-sm" aria-label={`About ${settings.title}`} />}
        >
          <Info aria-hidden />
        </HoverCardTrigger>
        <HoverCardContent side="bottom" align="end" className="grid w-72 gap-1.5">
          <span className="font-bold">{settings.title}</span>
          <span className="text-muted-foreground">{settings.description}</span>
          <code className="text-muted-foreground font-mono text-xs">treetime {command}</code>
        </HoverCardContent>
      </HoverCard>
    </div>
  );
}

function DetectedFacts({ command, facts }: { command: AppCommand; facts: readonly DetectedFact[] }) {
  if (facts.length === 0) {
    return null;
  }

  return (
    <dl className="grid grid-cols-[auto_minmax(0,1fr)_auto] items-baseline gap-x-3 gap-y-1 text-sm">
      {facts.map((fact) => (
        <DetectedFactRow key={fact.label} command={command} fact={fact} />
      ))}
    </dl>
  );
}

function DetectedFactRow({ command, fact }: { command: AppCommand; fact: DetectedFact }) {
  const spec = useMemo(
    () =>
      commandSettings(command).settings.find(
        (candidate) => candidate.key === fact.setting && candidate.role !== "input",
      ),
    [command, fact.setting],
  );

  return (
    <>
      <dt className="text-muted-foreground">{fact.label}</dt>
      <dd className="min-w-0 break-words">
        {fact.value}
        {fact.note !== undefined && <span className="text-muted-foreground ml-1.5 text-xs">({fact.note})</span>}
      </dd>
      <dd>
        {spec !== undefined && (
          <Popover>
            <PopoverTrigger
              render={
                <Button
                  type="button"
                  variant="link"
                  size="xs"
                  className="h-auto p-0"
                  aria-label={`Change ${fact.label}`}
                />
              }
            >
              Change
            </PopoverTrigger>
            <PopoverContent align="end" className="grid w-80 gap-2">
              <span className="font-bold">{spec.label}</span>
              <SettingControl command={command} spec={spec} label={spec.label} className="w-full" />
              <SettingHelp spec={spec} />
            </PopoverContent>
          </Popover>
        )}
      </dd>
    </>
  );
}

function AllSettings({ command, config }: { command: AppCommand; config: SparseConfig }) {
  const open = useDraftStore((state) => state.draft.view === "all");
  const update = useDraftStore((state) => state.update);
  const total = commandSettings(command).settings.filter((spec) => !spec.set_by_app).length;
  const onOpenChange = useCallback((next: boolean) => update({ view: next ? "all" : "main" }), [update]);
  const openAll = useCallback(() => update({ view: "all" }), [update]);

  return (
    <>
      <Button type="button" variant="outline" size="sm" className="justify-self-start" onClick={openAll}>
        <SlidersHorizontal aria-hidden />
        All settings
        <span className="text-muted-foreground">{total}</span>
      </Button>
      <Sheet open={open} onOpenChange={onOpenChange}>
        <SheetContent side="right" className="w-[min(100%,48rem)] gap-0 p-0 data-[side=right]:sm:max-w-3xl">
          <SheetHeader className="border-b">
            <SheetTitle>All settings</SheetTitle>
            <SheetDescription>{commandSettings(command).title}</SheetDescription>
          </SheetHeader>
          <div className="min-h-0 flex-1 overflow-y-auto overscroll-contain">
            <AllSettingsPanel command={command} config={config} />
          </div>
        </SheetContent>
      </Sheet>
    </>
  );
}

function Section({ title, hint, children }: { title: string; hint?: string; children: ReactNode }) {
  return (
    <section aria-label={title} className="grid gap-2.5">
      <div className="flex items-baseline gap-2 border-b pb-1.5">
        <h2 className="font-heading text-sm font-bold">{title}</h2>
        {hint !== undefined && <span className="text-muted-foreground text-xs">{hint}</span>}
      </div>
      {children}
    </section>
  );
}
