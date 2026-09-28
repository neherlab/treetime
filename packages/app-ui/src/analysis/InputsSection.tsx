import type { AppCommand, InputFacts, InputSlot } from "@neherlab/app-contracts";
import { useCallback, useState } from "react";
import { useWatch } from "react-hook-form";

import { formatBytes } from "../format";
import { useLocalFiles } from "../platform";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { baseName, pathList, slotFactsText, slotProblem } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { Button, cn } from "../ui";
import { ExamplesPanel } from "./ExamplesPanel";
import type { FormConfig } from "./formValues";
import { PathPicker, useFileDrop } from "./PathPicker";

export function InputsSection({ command, facts }: { command: AppCommand; facts: InputFacts | undefined }) {
  const [showExamples, setShowExamples] = useState(false);
  const localFiles = useLocalFiles();
  const toggleExamples = useCallback(() => setShowExamples((shown) => !shown), []);
  const closeExamples = useCallback(() => setShowExamples(false), []);

  return (
    <div className="grid grid-cols-1 gap-2">
      <div className="flex flex-wrap items-center gap-2">
        <span className="text-ink-faint">
          {localFiles === null
            ? "Files are uploaded into the run folder on the server"
            : "Files are read from this computer"}
        </span>
        <Button type="button" variant="outline" size="sm" className="ml-auto" onClick={toggleExamples}>
          {showExamples ? "Hide examples" : "Examples and earlier inputs"}
        </Button>
      </div>
      {showExamples && <ExamplesPanel command={command} close={closeExamples} />}
      {COMMAND_SETTINGS[command].inputs.map((input) => (
        <InputSlotRow key={input.kind} command={command} slot={input} facts={facts} />
      ))}
    </div>
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
  const drop = useFileDrop(command, slot.kind, slot.list);
  const paths = pathList(value);
  const factsText = slotFactsText(slot.kind, facts, COMMAND_SETTINGS[command].uses_dates);
  const problem = slotProblem(slot.kind, facts);

  return (
    <div
      onDragOver={drop.onDragOver}
      onDragLeave={drop.onDragLeave}
      onDrop={drop.onDrop}
      className={cn(
        "bg-surface-1 grid grid-cols-[7.5rem_minmax(0,1fr)_auto] items-center gap-3 rounded-lg border px-3 py-2.5",
        paths.length > 0 ? "border-line" : "border-line-strong border-dashed",
        drop.over && "border-accent bg-accent-subtle",
      )}
    >
      <div className="font-bold">
        {slot.label}
        <small className="text-ink-faint block text-xs font-normal">
          {slot.need === "required" ? slot.formats : `${slot.formats}, ${slot.need}`}
        </small>
      </div>
      {paths.length > 0 ? (
        <div className="min-w-0">
          <div className="truncate font-mono text-xs" title={paths.join("\n")}>
            {source?.label ?? paths.map(baseName).join(", ")}
            {source?.size !== null && source?.size !== undefined && (
              <span className="text-ink-faint ml-2">{formatBytes(source.size)}</span>
            )}
          </div>
          {problem === null ? (
            <div className="text-ink-muted text-xs">{factsText ?? "Checking..."}</div>
          ) : (
            <div className="text-signal-danger text-xs">{problem}</div>
          )}
        </div>
      ) : (
        <div className="text-ink-faint">Drop a file here or choose one</div>
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
