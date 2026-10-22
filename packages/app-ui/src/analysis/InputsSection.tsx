import type { AppCommand, InputFacts, InputSlot } from "@neherlab/app-contracts";
import { cn } from "cn";
import { useCallback, useMemo, useState } from "react";
import { useWatch } from "react-hook-form";
import Upload from "~icons/lucide/upload";

import { type BrowserBackend, useBrowserBackend } from "../browser-backend";
import { pressedChoice } from "../components/toggleChoice";
import { formatBytes } from "../format";
import type { Host } from "../host";
import { useHost } from "../host-context";
import { commandSettings } from "../settings/catalog";
import { baseName, pathList, slotFactsText, slotProblems } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { ToggleGroup, ToggleGroupItem } from "../ui/toggle-group";
import { EXAMPLES_PANEL_KINDS, EXAMPLES_PANELS, ExamplesPanel, type ExamplesPanelKind } from "./ExamplesPanel";
import type { FormConfig } from "./formValues";
import { PathPicker, useFileDrop } from "./PathPicker";

export function InputsSection({ command, facts }: { command: AppCommand; facts: InputFacts | undefined }) {
  const [panel, setPanel] = useState<ExamplesPanelKind | null>(null);
  const host = useHost();
  const browserBackend = useBrowserBackend();
  const pressed = useMemo(() => (panel === null ? [] : [panel]), [panel]);
  const closePanel = useCallback(() => setPanel(null), []);

  const onPanelChange = useCallback(
    (next: string[]) => setPanel(pressedChoice(EXAMPLES_PANEL_KINDS, next) ?? null),
    [],
  );

  return (
    <div className="grid grid-cols-1 gap-2">
      <div className="flex flex-wrap items-center gap-2">
        <span className="text-muted-foreground text-sm">{inputsNote(host, browserBackend)}</span>
        <ToggleGroup
          aria-label="Examples"
          variant="outline"
          size="sm"
          className="ml-auto"
          value={pressed}
          onValueChange={onPanelChange}
        >
          {EXAMPLES_PANEL_KINDS.map((kind) => {
            const { title, icon: Icon } = EXAMPLES_PANELS[kind];

            return (
              <ToggleGroupItem key={kind} value={kind}>
                <Icon aria-hidden />
                {title}
              </ToggleGroupItem>
            );
          })}
        </ToggleGroup>
      </div>
      {panel !== null && <ExamplesPanel kind={panel} command={command} close={closePanel} />}
      {commandSettings(command).inputs.map((input) => (
        <InputSlotRow key={input.kind} command={command} slot={input} facts={facts} />
      ))}
    </div>
  );
}

function inputsNote(host: Host | null, browserBackend: BrowserBackend | null): string {
  if (host !== null) {
    return "Files are read from this computer";
  }

  return browserBackend === null
    ? "Files are uploaded into the run folder on the server"
    : "Files and runs stay in this browser tab until it reloads or closes";
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
          "bg-card grid grid-cols-[7.5rem_minmax(0,1fr)_auto] items-center gap-3 rounded-lg border px-3 py-2.5 transition-colors",
          paths.length === 0 && "border-input border-dashed",
          isDragActive && "border-primary bg-accent",
        ),
      })}
    >
      <input {...getInputProps()} aria-label={`Dropped ${slot.label} files`} />
      <div className="grid gap-0.5">
        <span className="font-bold">{slot.label}</span>
        <span className="text-muted-foreground text-xs">
          {slot.need === "required" ? slot.formats : `${slot.formats}, ${slot.need}`}
        </span>
      </div>
      {paths.length > 0 ? (
        <div className="grid min-w-0 gap-0.5">
          <div className="truncate font-mono text-xs" title={paths.join("\n")}>
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
