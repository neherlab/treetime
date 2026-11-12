import type { AppCommand, SparseConfig } from "@neherlab/app-contracts";
import { Badge } from "@neherlab/design/ui/badge";
import { Checkbox } from "@neherlab/design/ui/checkbox";
import { Empty, EmptyDescription } from "@neherlab/design/ui/empty";
import { InputGroup, InputGroupAddon, InputGroupInput } from "@neherlab/design/ui/input-group";
import { Label } from "@neherlab/design/ui/label";
import { useCallback } from "react";
import Search from "~icons/lucide/search";

import { commandSettings } from "../settings/catalog";
import { isChanged } from "../settings/config";
import { settingSections, type SettingSectionView } from "../settings/groups";
import { matchingSpecs } from "../settings/search";
import { useDraftStore } from "../store/draft";
import { SettingField } from "./SettingField";

export function AllSettingsPanel({ command, config }: { command: AppCommand; config: SparseConfig }) {
  const search = useDraftStore((state) => state.draft.search);
  const changedOnly = useDraftStore((state) => state.draft.changed_only);
  const update = useDraftStore((state) => state.update);

  const onSearch = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => update({ search: event.target.value }),
    [update],
  );

  const onChangedOnly = useCallback((checked: boolean) => update({ changed_only: checked }), [update]);

  return (
    <div className="@container">
      <div className="flex flex-wrap items-center gap-2.5 border-b px-3.5 py-2.5">
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
      </div>
      <AllSettings command={command} config={config} search={search} changedOnly={changedOnly} />
    </div>
  );
}

function AllSettings({
  command,
  config,
  search,
  changedOnly,
}: {
  command: AppCommand;
  config: SparseConfig;
  search: string;
  changedOnly: boolean;
}) {
  const settings = commandSettings(command);
  const sections = settingSections(settings, matchingSpecs(settings.settings, config, search, changedOnly));
  const allSections = settingSections(settings, settings.settings);

  if (sections.length === 0) {
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
      <nav aria-label="Setting sections" className="flex flex-wrap gap-1 border-b px-3.5 py-2">
        {sections.map(({ purpose, title }) => {
          const changed = changedCount(
            config,
            allSections.find((candidate) => candidate.purpose === purpose),
          );

          return (
            <Badge
              key={purpose}
              variant="outline"
              render={
                <a href={`#section-${purpose}`} aria-label={changed > 0 ? `${title}, ${changed} changed` : title} />
              }
            >
              {title}
              {changed > 0 && <span className="text-primary">{changed} changed</span>}
            </Badge>
          );
        })}
      </nav>
      {sections.map((section) => (
        <SettingSection key={section.purpose} command={command} config={config} section={section} />
      ))}
    </div>
  );
}

function SettingSection({
  command,
  config,
  section,
}: {
  command: AppCommand;
  config: SparseConfig;
  section: SettingSectionView;
}) {
  const specs = section.groups.flatMap((group) => group.specs);
  const changed = changedCount(config, section);
  const headed = section.groups.length > 1 || section.groups.some((group) => group.group !== section.title);

  return (
    <section
      id={`section-${section.purpose}`}
      aria-label={section.title}
      className="scroll-mt-4 border-b pb-2 last:border-b-0"
    >
      <h3 className="flex items-baseline gap-2 px-3.5 pt-3 pb-1.5 font-bold">
        {section.title}
        <span className="text-muted-foreground text-xs font-normal">
          {specs.length} {specs.length === 1 ? "setting" : "settings"}
        </span>
        {changed > 0 && <span className="text-primary text-xs font-normal">{changed} changed</span>}
      </h3>
      {section.groups.map(({ group, specs: members }) => (
        <div key={group}>
          {headed && (
            <h4 className="text-muted-foreground px-3.5 pt-2 pb-1 text-xs font-bold tracking-wide uppercase">
              {group}
            </h4>
          )}
          {members.map((spec) => (
            <SettingField key={spec.key} command={command} spec={spec} config={config} />
          ))}
        </div>
      ))}
    </section>
  );
}

function changedCount(config: SparseConfig, section: SettingSectionView | undefined): number {
  return (section?.groups ?? []).flatMap((group) => group.specs).filter((spec) => isChanged(config, spec)).length;
}
