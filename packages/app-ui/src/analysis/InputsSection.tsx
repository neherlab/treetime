import type { AppCommand, InputFacts, InputSlot } from "@neherlab/app-contracts";
import { BookOpen, Upload } from "lucide-react";
import { useCallback, useState } from "react";
import { useWatch } from "react-hook-form";

import { formatBytes } from "../format";
import { useLocalFiles } from "../platform";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { baseName, pathList, slotFactsText, slotProblem } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { Button } from "../ui/button";
import { cn } from "../ui/cn";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "../ui/collapsible";
import { ExamplesPanel } from "./ExamplesPanel";
import type { FormConfig } from "./formValues";
import { PathPicker, useFileDrop } from "./PathPicker";

export function InputsSection({ command, facts }: { command: AppCommand; facts: InputFacts | undefined }) {
  const [showExamples, setShowExamples] = useState(false);
  const localFiles = useLocalFiles();
  const closeExamples = useCallback(() => setShowExamples(false), []);

  return (
    <Collapsible open={showExamples} onOpenChange={setShowExamples} className="grid grid-cols-1 gap-2">
      <div className="flex flex-wrap items-center gap-2">
        <span className="text-muted-foreground text-sm">
          {localFiles === null
            ? "Files are uploaded into the run folder on the server"
            : "Files are read from this computer"}
        </span>
        <CollapsibleTrigger render={<Button type="button" variant="outline" size="sm" className="ml-auto" />}>
          <BookOpen aria-hidden />
          {showExamples ? "Hide examples" : "Examples and earlier inputs"}
        </CollapsibleTrigger>
      </div>
      <CollapsibleContent>
        <ExamplesPanel command={command} close={closeExamples} />
      </CollapsibleContent>
      {COMMAND_SETTINGS[command].inputs.map((input) => (
        <InputSlotRow key={input.kind} command={command} slot={input} facts={facts} />
      ))}
    </Collapsible>
  );
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
  const value = useWatch<FormConfig>({ name: slot.kind });
  const source = useDraftStore((state) => state.sources[slot.kind]);
  const { getRootProps, getInputProps, isDragActive } = useFileDrop(command, slot.kind, slot.list);
  const paths = pathList(value);
  const factsText = slotFactsText(slot.kind, facts, COMMAND_SETTINGS[command].uses_dates);
  const problem = slotProblem(slot.kind, facts);

  return (
    <div
      {...getRootProps({
        className: cn(
          "bg-card grid grid-cols-[7.5rem_minmax(0,1fr)_auto] items-center gap-3 rounded-lg border px-3 py-2.5 transition-colors",
          paths.length === 0 && "border-input border-dashed",
          isDragActive && "border-primary bg-accent",
        ),
      })}
    >
      <input {...getInputProps()} />
      <div className="grid gap-0.5">
        <span className="font-medium">{slot.label}</span>
        <span className="text-muted-foreground text-xs">
          {slot.need === "required" ? slot.formats : `${slot.formats}, ${slot.need}`}
        </span>
      </div>
      {paths.length > 0 ? (
        <div className="grid min-w-0 gap-0.5">
          <div className="truncate font-mono text-xs" title={paths.join("\n")}>
            {source?.label ?? paths.map(baseName).join(", ")}
            {source?.size !== null && source?.size !== undefined && (
              <span className="text-muted-foreground ml-2">{formatBytes(source.size)}</span>
            )}
          </div>
          {problem === null ? (
            <div className="text-muted-foreground text-xs">{factsText ?? "Checking..."}</div>
          ) : (
            <div className="text-destructive text-xs">{problem}</div>
          )}
        </div>
      ) : (
        <div className="text-muted-foreground flex items-center gap-2 text-sm">
          <Upload aria-hidden className="size-4" />
          Drop a file here or choose one
        </div>
      )}
      <PathPicker
        command={command}
        settingKey={slot.kind}
        title={slot.label}
        extensions={slot.extensions}
        list={slot.list}
        filled={paths.length > 0}
      />
    </div>
  );
}
