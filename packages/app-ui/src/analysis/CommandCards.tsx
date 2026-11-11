import type {
  AppCommand,
  CommandSettings,
  InputKind,
  InputNeed,
  InputSlot,
  SparseConfig,
} from "@neherlab/app-contracts";
import { Field, FieldContent, FieldDescription, FieldLabel, FieldTitle } from "@neherlab/design/ui/field";
import { RadioGroup, RadioGroupItem } from "@neherlab/design/ui/radio-group";
import { cn } from "cn";
import { useCallback, type ComponentType, type SVGProps } from "react";
import Calendar from "~icons/lucide/calendar";
import Dna from "~icons/lucide/dna";
import TreeDeciduous from "~icons/lucide/tree-deciduous";

import { COMMANDS } from "../settings/catalog";
import { carryOverDraft } from "../settings/commands";
import { useDraftStore } from "../store/draft";

export function CommandCards({ command, config }: { command: AppCommand; config: SparseConfig }) {
  const select = useCallback(
    (value: string) => {
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

  return (
    <FieldLabel
      htmlFor={id}
      className="bg-card has-data-checked:border-primary has-data-checked:ring-primary dark:has-data-checked:border-primary has-data-checked:ring-1"
    >
      <Field orientation="horizontal" className="h-full items-start">
        <FieldContent className="h-full">
          <FieldTitle className="w-full justify-between">
            {settings.title}
            <code className="text-muted-foreground font-mono text-xs font-normal">{command}</code>
          </FieldTitle>
          <FieldDescription>{settings.description}</FieldDescription>
          <InputBadges inputs={settings.inputs} />
        </FieldContent>
        <RadioGroupItem value={command} id={id} className="sr-only" />
      </Field>
    </FieldLabel>
  );
}

const INPUT_ICONS: Record<InputKind, ComponentType<SVGProps<SVGSVGElement>>> = {
  tree: TreeDeciduous,
  metadata: Calendar,
  alignment: Dna,
};

const INPUT_ORDER: readonly InputKind[] = ["tree", "metadata", "alignment"];

const NEED_TEXT: Record<InputNeed, string> = {
  required: "Required",
  recommended: "Recommended",
  optional: "Optional",
};

function InputBadges({ inputs }: { inputs: readonly InputSlot[] }) {
  const ordered = INPUT_ORDER.flatMap((kind) => inputs.filter((input) => input.kind === kind));

  return (
    <ul aria-label="Inputs" className="mt-auto flex gap-1.5 pt-2">
      {ordered.map((input) => {
        const Icon = INPUT_ICONS[input.kind];
        const required = input.need === "required";
        const text = `${NEED_TEXT[input.need]} input: ${input.label.toLowerCase()}`;

        return (
          <li
            key={input.kind}
            title={text}
            className={cn(
              "flex items-center rounded-sm border p-1",
              required
                ? "border-foreground/30 text-foreground"
                : "border-muted-foreground/40 text-muted-foreground border-dashed",
            )}
          >
            <Icon aria-hidden className="size-4" />
            <span className="sr-only">{text}</span>
          </li>
        );
      })}
    </ul>
  );
}
