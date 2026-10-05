import type { AppCommand, CommandSettings, SparseConfig } from "@neherlab/app-contracts";
import { useCallback } from "react";

import { COMMANDS, commandSettings } from "../settings/catalog";
import { carryOverConfig } from "../settings/config";
import { useDraftStore } from "../store/draft";
import { Field, FieldContent, FieldDescription, FieldLabel, FieldTitle } from "../ui/field";
import { RadioGroup, RadioGroupItem } from "../ui/radio-group";

export function CommandCards({ command, config }: { command: AppCommand; config: SparseConfig }) {
  const select = useCallback(
    (value: string) => {
      const target = COMMANDS.find((candidate) => candidate.command === value)?.command;

      if (target === undefined || target === command) {
        return;
      }

      const targetSpecs = commandSettings(target).settings;
      const keys = new Set(targetSpecs.map((spec) => spec.key));
      const { draft, load } = useDraftStore.getState();

      load({
        ...draft,
        command: target,
        config: carryOverConfig(targetSpecs, commandSettings(command).settings, config),
        sources: Object.fromEntries(Object.entries(draft.sources).filter(([key]) => keys.has(key))),
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
      {COMMANDS.map((candidate) => (
        <CommandCard key={candidate.command} settings={candidate} />
      ))}
    </RadioGroup>
  );
}

function CommandCard({ settings }: { settings: CommandSettings }) {
  const { command } = settings;
  const id = `command-${command}`;

  const needs = settings.inputs
    .flatMap((input) => (input.need === "required" ? [input.label.toLowerCase()] : []))
    .join(" and ");

  return (
    <FieldLabel htmlFor={id} className="bg-card">
      <Field orientation="horizontal" className="items-start">
        <FieldContent>
          <FieldTitle className="w-full justify-between">
            {settings.title}
            <code className="text-muted-foreground font-mono text-xs font-normal">{command}</code>
          </FieldTitle>
          <FieldDescription>{settings.description}</FieldDescription>
          <FieldDescription className="text-xs">Needs {needs}</FieldDescription>
        </FieldContent>
        <RadioGroupItem value={command} id={id} />
      </Field>
    </FieldLabel>
  );
}
