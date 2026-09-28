import type { AppCommand } from "@neherlab/app-contracts";
import { useCallback } from "react";

import { COMMAND_SETTINGS } from "../settings/catalog";
import { APP_COMMANDS, COMMAND_INFO } from "../settings/commands";
import { carryOverConfig } from "../settings/config";
import type { JsonObject } from "../settings/json";
import { useDraftStore } from "../store/draft";
import { Field, FieldContent, FieldDescription, FieldLabel, FieldTitle } from "../ui/field";
import { RadioGroup, RadioGroupItem } from "../ui/radio-group";

export function CommandCards({ command, config }: { command: AppCommand; config: JsonObject }) {
  const select = useCallback(
    (value: string) => {
      const target = APP_COMMANDS.find((candidate) => candidate === value);

      if (target === undefined || target === command) {
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
    <RadioGroup
      aria-label="Analysis"
      value={command}
      onValueChange={select}
      className="grid-cols-2 gap-2 @2xl:grid-cols-3"
    >
      {APP_COMMANDS.map((candidate) => (
        <CommandCard key={candidate} command={candidate} />
      ))}
    </RadioGroup>
  );
}

function CommandCard({ command }: { command: AppCommand }) {
  const info = COMMAND_INFO[command];
  const id = `command-${command}`;

  const needs = COMMAND_SETTINGS[command].inputs
    .flatMap((input) => (input.need === "required" ? [input.label.toLowerCase()] : []))
    .join(" and ");

  return (
    <FieldLabel htmlFor={id} className="bg-card">
      <Field orientation="horizontal" className="items-start">
        <FieldContent>
          <FieldTitle className="w-full justify-between">
            {info.label}
            <code className="text-muted-foreground font-mono text-xs font-normal">{command}</code>
          </FieldTitle>
          <FieldDescription>{info.description}</FieldDescription>
          <FieldDescription className="text-xs">Needs {needs}</FieldDescription>
        </FieldContent>
        <RadioGroupItem value={command} id={id} />
      </Field>
    </FieldLabel>
  );
}
