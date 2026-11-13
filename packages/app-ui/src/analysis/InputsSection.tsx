import type { AppCommand, InputFacts, InputSlot } from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { cn } from "cn";
import { useCallback, useState } from "react";
import { useWatch } from "react-hook-form";
import History from "~icons/lucide/history";
import Upload from "~icons/lucide/upload";

import { formatBytes } from "../format";
import type { Host } from "../host";
import { useHost } from "../host-context";
import { commandSettings } from "../settings/catalog";
import { baseName, pathList, slotFactsText, slotProblems } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { EarlierInputsPanel } from "./EarlierInputsPanel";
import type { FormConfig } from "./formValues";
import { PathPicker, useFileDrop } from "./PathPicker";

export function InputsSection({ command, facts }: { command: AppCommand; facts: InputFacts | undefined }) {
  const [earlierOpen, setEarlierOpen] = useState(false);
  const host = useHost();
  const closeEarlier = useCallback(() => setEarlierOpen(false), []);
  const toggleEarlier = useCallback(() => setEarlierOpen((open) => !open), []);

  return (
    <div className="grid grid-cols-1 gap-2">
      <div className="flex flex-wrap items-center gap-2">
        <span className="text-muted-foreground text-sm">{inputsNote(host)}</span>
        <Button
          type="button"
          variant="outline"
          size="sm"
          className="ml-auto"
          aria-pressed={earlierOpen}
          onClick={toggleEarlier}
        >
          <History aria-hidden />
          Inputs of earlier runs
        </Button>
      </div>
      {earlierOpen && <EarlierInputsPanel command={command} close={closeEarlier} />}
      {commandSettings(command).inputs.map((input) => (
        <InputSlotRow key={input.kind} command={command} slot={input} facts={facts} />
      ))}
    </div>
  );
}

function inputsNote(host: Host | null): string {
  return host === null
    ? "Files and runs stay in this browser tab. A reload, a closed tab, or a crash of the computation loses them."
    : "Files are read from this computer";
}

function InputSlotRow({
  command,
  slot,
  facts,
}: {
  command: AppCommand;
  slot: InputSlot;
  facts: InputFacts | undefined;
}) {
  const watched = useWatch<FormConfig>({ name: slot.kind });
  const source = useDraftStore((state) => state.draft.sources[slot.kind]);
  const { getRootProps, getInputProps, isDragActive } = useFileDrop(command, slot.kind, slot.list);
  const paths = pathList(watched);
  const factsText = slotFactsText(slot.kind, facts, commandSettings(command).uses_dates);
  const problems = slotProblems(slot.kind, facts);

  return (
    <div
      {...getRootProps({
        role: "group",
        "aria-label": slot.label,
        className: cn(
          "bg-card grid grid-cols-[minmax(0,1fr)_auto] items-center gap-x-3 gap-y-1 rounded-lg border px-3 py-2.5 transition-colors",
          paths.length === 0 && "border-input border-dashed",
          isDragActive && "border-primary bg-accent",
        ),
      })}
    >
      <input {...getInputProps({ className: "absolute" })} aria-label={`Dropped ${slot.label} files`} />
      <div className="flex min-w-0 flex-wrap items-baseline gap-x-2">
        <span className="font-bold">{slot.label}</span>
        <span className="text-muted-foreground text-xs">
          {slot.need === "required" ? slot.formats : `${slot.formats}, ${slot.need}`}
        </span>
      </div>
      <div className="row-span-2">
        <PathPicker
          command={command}
          settingKey={slot.kind}
          title={slot.label}
          extensions={slot.extensions}
          list={slot.list}
          filled={paths.length > 0}
        />
      </div>
      {paths.length > 0 ? (
        <div className="grid min-w-0 gap-0.5">
          <div className="font-mono text-xs break-all" title={paths.join("\n")}>
            {source?.label ?? paths.map(baseName).join(", ")}
            {source?.size !== undefined && (
              <span className="text-muted-foreground ml-2">{formatBytes(source.size)}</span>
            )}
          </div>
          {(factsText ?? (problems.length === 0 ? ["Checking..."] : [])).map((line) => (
            <div key={line} className="text-muted-foreground text-xs">
              {line}
            </div>
          ))}
          {problems.map((problem) => (
            <div key={problem} className="text-destructive text-xs">
              {problem}
            </div>
          ))}
        </div>
      ) : (
        <div className="text-muted-foreground flex items-center gap-2 text-sm">
          <Upload aria-hidden className="size-4" />
          {slot.list ? "Drop files here or choose them" : "Drop a file here or choose one"}
        </div>
      )}
    </div>
  );
}
