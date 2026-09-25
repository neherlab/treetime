import type { AppCommand, InputFactsResult } from "@neherlab/app-contracts";
import { useCallback, useMemo } from "react";

import { COMMAND_SETTINGS } from "../settings/catalog";
import { isChanged } from "../settings/config";
import { groupedSpecs } from "../settings/groups";
import type { JsonObject } from "../settings/json";
import { matchingSpecs } from "../settings/search";
import { useDraftStore, type SettingsView } from "../store/draft";
import { Segmented } from "../ui";
import { MainSettings } from "./MainSettings";
import { SettingField } from "./SettingField";

export function SettingsPanel({
  command,
  config,
  facts,
}: {
  command: AppCommand;
  config: JsonObject;
  facts: InputFactsResult | undefined;
}) {
  const view = useDraftStore((state) => state.view);
  const search = useDraftStore((state) => state.search);
  const changedOnly = useDraftStore((state) => state.changedOnly);
  const update = useDraftStore((state) => state.update);
  const total = COMMAND_SETTINGS[command].specs.filter((spec) => spec.pathRole !== "output").length;

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

  const onChangedOnly = useCallback(
    (event: React.ChangeEvent<HTMLInputElement>) => update({ changedOnly: event.target.checked }),
    [update],
  );

  return (
    <div className="border-line bg-surface-1 rounded-lg border">
      <div className="border-line flex flex-wrap items-center gap-2.5 border-b px-3.5 py-2.5">
        <Segmented label="Settings view" value={view} onChange={onView} options={views} />
        {view === "all" ? (
          <>
            <input
              type="search"
              value={search}
              onChange={onSearch}
              placeholder="Search settings by name, flag or description"
              aria-label="Search settings"
              className="border-line bg-surface-2 min-w-44 flex-1 rounded-md border px-2.5 py-1"
            />
            <label className="text-ink-muted inline-flex items-center gap-1.5">
              <input type="checkbox" checked={changedOnly} onChange={onChangedOnly} />
              Changed only
            </label>
          </>
        ) : (
          <span className="text-ink-faint text-xs">
            Every setting is also under All settings, with its CLI flag and help text.
          </span>
        )}
      </div>
      {view === "main" ? (
        <MainSettings command={command} config={config} facts={facts} />
      ) : (
        <AllSettings command={command} config={config} search={search} changedOnly={changedOnly} />
      )}
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
  config: JsonObject;
  search: string;
  changedOnly: boolean;
}) {
  const specs = COMMAND_SETTINGS[command].specs;
  const groups = groupedSpecs(matchingSpecs(specs, config, search, changedOnly));
  const allGroups = groupedSpecs(specs);
  const filtering = search.trim() !== "" || changedOnly;

  if (groups.length === 0) {
    return (
      <p className="text-ink-muted p-10 text-center">
        {changedOnly && search.trim() === ""
          ? "No setting differs from its default."
          : `No setting matches "${search}".`}
      </p>
    );
  }

  return (
    <div>
      <nav aria-label="Setting groups" className="border-line flex flex-wrap gap-1 border-b px-3.5 py-2">
        {groups.map(([group]) => {
          const changed =
            allGroups.find(([candidate]) => candidate === group)?.[1].filter((spec) => isChanged(config, spec))
              .length ?? 0;

          return (
            <a
              key={group}
              href={`#group-${groupAnchor(group)}`}
              className="border-line text-ink-muted hover:border-line-strong rounded-full border px-2 py-0.5 text-xs"
            >
              {group}
              {changed > 0 && ` - ${changed} changed`}
            </a>
          );
        })}
      </nav>
      {groups.map(([group, members]) => {
        const editable = members.filter((spec) => spec.pathRole !== "output");
        const outputs = members.filter((spec) => spec.pathRole === "output");
        const changed = members.filter((spec) => isChanged(config, spec)).length;

        return (
          <details
            key={group}
            id={`group-${groupAnchor(group)}`}
            open={filtering || group !== "Outputs"}
            className="border-line border-b last:border-b-0"
          >
            <summary className="flex cursor-pointer items-center gap-2 px-3.5 py-2.5 font-bold">
              {group}
              <span className="text-ink-faint text-xs font-normal">
                {members.length} {members.length === 1 ? "setting" : "settings"}
              </span>
              {changed > 0 && <span className="text-accent text-xs">{changed} changed</span>}
            </summary>
            {editable.map((spec) => (
              <SettingField key={spec.key} command={command} spec={spec} config={config} />
            ))}
            {outputs.length > 0 && (
              <details className="px-7.5 pt-1.5 pb-2.5">
                <summary className="text-ink-muted cursor-pointer text-[0.8125rem]">
                  {outputs.length} file paths, set by the app for each run
                </summary>
                {outputs.map((spec) => (
                  <SettingField key={spec.key} command={command} spec={spec} config={config} />
                ))}
              </details>
            )}
          </details>
        );
      })}
    </div>
  );
}

function groupAnchor(group: string): string {
  return group.toLowerCase().replaceAll(/\W+/gu, "-");
}
