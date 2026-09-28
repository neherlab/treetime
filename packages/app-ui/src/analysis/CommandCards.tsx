import type { AppCommand } from "@neherlab/app-contracts";
import { useCallback } from "react";

import { COMMAND_SETTINGS } from "../settings/catalog";
import { APP_COMMANDS, COMMAND_INFO } from "../settings/commands";
import { carryOverConfig } from "../settings/config";
import type { JsonObject } from "../settings/json";
import { useDraftStore } from "../store/draft";

export function CommandCards({ command, config }: { command: AppCommand; config: JsonObject }) {
  const select = useCallback(
    (target: AppCommand) => {
      if (target === command) {
        return;
      }

      const targetSpecs = COMMAND_SETTINGS[target].specs;
      const keys = new Set(targetSpecs.map((spec) => spec.key));
      const { sources, load } = useDraftStore.getState();

      load({
        command: target,
        config: carryOverConfig(targetSpecs, COMMAND_SETTINGS[command].specs, config),
        sources: Object.fromEntries(Object.entries(sources).filter(([key]) => keys.has(key))),
      });
    },
    [command, config],
  );

  return (
    <fieldset className="m-0 grid grid-cols-2 gap-2 border-0 p-0 @2xl:grid-cols-3">
      <legend className="sr-only">Analysis</legend>
      {APP_COMMANDS.map((candidate) => (
        <CommandCard key={candidate} command={candidate} pressed={candidate === command} select={select} />
      ))}
    </fieldset>
  );
}

function CommandCard({
  command,
  pressed,
  select,
}: {
  command: AppCommand;
  pressed: boolean;
  select: (command: AppCommand) => void;
}) {
  const info = COMMAND_INFO[command];
  const onClick = useCallback(() => select(command), [command, select]);

  const needs = COMMAND_SETTINGS[command].inputs
    .flatMap((input) => (input.need === "required" ? [input.label.toLowerCase()] : []))
    .join(" and ");

  return (
    <button
      type="button"
      aria-pressed={pressed}
      onClick={onClick}
      className="border-line bg-surface-1 hover:border-line-strong aria-pressed:border-accent aria-pressed:bg-accent-subtle aria-pressed:ring-accent grid gap-1 rounded-lg border px-3 py-2.5 text-left aria-pressed:ring-1 aria-pressed:ring-inset"
    >
      <span className="flex justify-between text-[0.9375rem] font-bold">
        {info.label}
        <code className="text-ink-faint font-mono text-xs font-normal">{command}</code>
      </span>
      <span className="text-ink-muted">{info.description}</span>
      <span className="text-ink-faint text-xs">Needs {needs}</span>
    </button>
  );
}
