import type { AppCommand, InputFacts, RunCheck } from "@neherlab/app-contracts";
import { useCallback, useMemo } from "react";
import ChevronRight from "~icons/lucide/chevron-right";
import Search from "~icons/lucide/search";

import { OptionToggle } from "../components/OptionToggle";
import { COMMAND_SETTINGS, groupedSpecs, type SettingSpec } from "../settings/catalog";
import { isChanged } from "../settings/config";
import type { JsonObject } from "../settings/json";
import { matchingSpecs } from "../settings/search";
import { useDraftStore } from "../store/draft";
import type { SettingsView } from "../store/draftSchema";
import { Badge } from "../ui/badge";
import { Card } from "../ui/card";
import { Checkbox } from "../ui/checkbox";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "../ui/collapsible";
import { Empty, EmptyDescription } from "../ui/empty";
import { InputGroup, InputGroupAddon, InputGroupInput } from "../ui/input-group";
import { Label } from "../ui/label";
import { MainSettings } from "./MainSettings";
import { SettingField } from "./SettingField";

export function SettingsPanel({
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
  const view = useDraftStore((state) => state.view);
  const search = useDraftStore((state) => state.search);
  const changedOnly = useDraftStore((state) => state.changedOnly);
  const update = useDraftStore((state) => state.update);
  const total = COMMAND_SETTINGS[command].specs.filter((spec) => spec.role !== "output").length;

  const views = useMemo(
    () => [
      { value: "main" as const, label: "Main settings" },
      { value: "all" as const, label: `All ${total} settings` },
    ],
    [total],
  );

  const onView = useCallback((next: SettingsView) => update({ view: next }), [update]);

  const onSearch = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => update({ search: event.target.value }),
    [update],
  );

  const onChangedOnly = useCallback((checked: boolean) => update({ changedOnly: checked }), [update]);

  return (
    <Card size="sm" className="gap-0 py-0">
      <div className="flex flex-wrap items-center gap-2.5 border-b px-3.5 py-2.5">
        <OptionToggle label="Settings view" value={view} onChange={onView} options={views} />
        {view === "all" ? (
          <>
            <InputGroup className="h-8 min-w-44 flex-1">
              <InputGroupAddon>
                <Search aria-hidden />
              </InputGroupAddon>
              <InputGroupInput
                type="search"
                value={search}
                onChange={onSearch}
                placeholder="Search settings by name, flag or description"
                aria-label="Search settings"
              />
            </InputGroup>
            <div className="flex items-center gap-2">
              <Checkbox id="settings-changed-only" checked={changedOnly} onCheckedChange={onChangedOnly} />
              <Label htmlFor="settings-changed-only">Changed only</Label>
            </div>
          </>
        ) : (
          <span className="text-muted-foreground text-xs">
            Every setting is also under All settings, with its CLI flag and help text.
          </span>
        )}
      </div>
      {view === "main" ? (
        <MainSettings command={command} config={config} facts={facts} checks={checks} />
      ) : (
        <AllSettings command={command} config={config} search={search} changedOnly={changedOnly} />
      )}
    </Card>
  );
}

function AllSettings({
  command,
  config,
  search,
  changedOnly,
}: {
  command: AppCommand;
  config: JsonObject;
  search: string;
  changedOnly: boolean;
}) {
  const settings = COMMAND_SETTINGS[command];
  const specs = settings.specs;
  const groups = groupedSpecs(settings, matchingSpecs(specs, config, search, changedOnly));
  const allGroups = groupedSpecs(settings, specs);
  const filtering = search.trim() !== "" || changedOnly;

  if (groups.length === 0) {
    return (
      <Empty className="py-10">
        <EmptyDescription>
          {changedOnly && search.trim() === ""
            ? "No setting differs from its default."
            : `No setting matches "${search}".`}
        </EmptyDescription>
      </Empty>
    );
  }

  return (
    <div>
      <nav aria-label="Setting groups" className="flex flex-wrap gap-1 border-b px-3.5 py-2">
        {groups.map(([group]) => {
          const changed =
            allGroups.find(([candidate]) => candidate === group)?.[1].filter((spec) => isChanged(config, spec))
              .length ?? 0;

          return (
            <Badge
              key={group}
              variant="outline"
              render={
                <a
                  href={`#group-${groupAnchor(group)}`}
                  aria-label={changed > 0 ? `${group}, ${changed} changed` : group}
                />
              }
            >
              {group}
              {changed > 0 && <span className="text-primary">{changed} changed</span>}
            </Badge>
          );
        })}
      </nav>
      {groups.map(([group, members]) => (
        <SettingGroup
          key={`${group}:${String(filtering)}`}
          command={command}
          config={config}
          group={group}
          members={members}
          filtering={filtering}
        />
      ))}
    </div>
  );
}

function SettingGroup({
  command,
  config,
  group,
  members,
  filtering,
}: {
  command: AppCommand;
  config: JsonObject;
  group: string;
  members: readonly SettingSpec[];
  filtering: boolean;
}) {
  const editable = members.filter((spec) => spec.role !== "output");
  const outputs = members.filter((spec) => spec.role === "output");
  const changed = members.filter((spec) => isChanged(config, spec)).length;

  return (
    <Collapsible
      id={`group-${groupAnchor(group)}`}
      defaultOpen={filtering || outputs.length === 0}
      className="scroll-mt-4 border-b last:border-b-0"
    >
      <CollapsibleTrigger className="group/trigger hover:bg-muted/50 flex w-full items-center gap-2 px-3.5 py-2.5 text-left font-bold">
        <ChevronRight
          aria-hidden
          className="text-muted-foreground size-4 transition-transform group-data-panel-open/trigger:rotate-90"
        />
        {group}
        <span className="text-muted-foreground text-xs font-normal">
          {members.length} {members.length === 1 ? "setting" : "settings"}
        </span>
        {changed > 0 && <span className="text-primary text-xs">{changed} changed</span>}
      </CollapsibleTrigger>
      <CollapsibleContent>
        {editable.map((spec) => (
          <SettingField key={spec.key} command={command} spec={spec} config={config} />
        ))}
        {outputs.length > 0 && (
          <Collapsible className="px-7.5 pt-1.5 pb-2.5">
            <CollapsibleTrigger className="text-muted-foreground hover:text-foreground text-xs">
              {outputs.length} file paths, set by the app for each run
            </CollapsibleTrigger>
            <CollapsibleContent>
              {outputs.map((spec) => (
                <SettingField key={spec.key} command={command} spec={spec} config={config} />
              ))}
            </CollapsibleContent>
          </Collapsible>
        )}
      </CollapsibleContent>
    </Collapsible>
  );
}

function groupAnchor(group: string): string {
  return group.toLowerCase().replaceAll(/\W+/gu, "-");
}
