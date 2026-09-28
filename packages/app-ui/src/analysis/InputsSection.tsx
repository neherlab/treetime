import type { AppCommand, InputFacts, InputSlot } from "@neherlab/app-contracts";
import { Upload } from "lucide-react";
import { useCallback, useMemo, useState } from "react";
import { useWatch } from "react-hook-form";

import { pressedChoice } from "../components/toggleChoice";
import { formatBytes } from "../format";
import { useLocalFiles } from "../platform";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { baseName, pathList, slotFactsText, slotProblem } from "../settings/inputs";
import { useDraftStore } from "../store/draft";
import { cn } from "../ui/cn";
import { ToggleGroup, ToggleGroupItem } from "../ui/toggle-group";
import { EXAMPLES_PANEL_KINDS, EXAMPLES_PANELS, ExamplesPanel, type ExamplesPanelKind } from "./ExamplesPanel";
import type { FormConfig } from "./formValues";
import { PathPicker, useFileDrop } from "./PathPicker";

export function InputsSection({ command, facts }: { command: AppCommand; facts: InputFacts | undefined }) {
  const [panel, setPanel] = useState<ExamplesPanelKind | null>(null);
  const localFiles = useLocalFiles();
  const pressed = useMemo(() => (panel === null ? [] : [panel]), [panel]);
  const closePanel = useCallback(() => setPanel(null), []);

  const onPanelChange = useCallback(
    (next: string[]) => setPanel(pressedChoice(EXAMPLES_PANEL_KINDS, next) ?? null),
    [],
  );

  return (
    <div className="grid grid-cols-1 gap-2">
      <div className="flex flex-wrap items-center gap-2">
        <span className="text-muted-foreground text-sm">
          {localFiles === null
            ? "Files are uploaded into the run folder on the server"
            : "Files are read from this computer"}
        </span>
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
        <span className="font-bold">{slot.label}</span>
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
