import { errorMessage } from "@neherlab/app-contracts";
import type { AppCommand, CodeLine, ConfigCode } from "@neherlab/app-contracts";
import { useCallback, useMemo, useState } from "react";
import FileUp from "~icons/lucide/file-up";

import { CopyButton } from "../components/CopyButton";
import { OptionToggle } from "../components/OptionToggle";
import { Panel } from "../components/Panel";
import { commandSwitchNote } from "../settings/commands";
import { useDraftStore } from "../store/draft";
import type { CodeFormat } from "../store/draftSchema";
import { Button } from "../ui/button";
import { cn } from "../ui/cn";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "../ui/collapsible";
import { Spinner } from "../ui/spinner";
import { Textarea } from "../ui/textarea";
import { useToastManager } from "../ui/toast";
import { useConfigLoader } from "./useConfigLoader";

const CODE_FORMATS: ReadonlyArray<{ value: CodeFormat; label: string }> = [
  { value: "cli", label: "CLI" },
  { value: "yaml", label: "YAML" },
];

export function CodePanel({ command, code }: { command: AppCommand; code: ConfigCode | null }) {
  const format = useDraftStore((state) => state.codeFormat);
  const update = useDraftStore((state) => state.update);
  const [importing, setImporting] = useState(false);

  const lines = useMemo(() => (code === null ? [] : format === "cli" ? code.command_line : code.yaml), [code, format]);
  const text = code === null ? "" : format === "cli" ? code.command_line_text : code.yaml_text;
  const keyed = useMemo(() => keyedLines(lines), [lines]);

  const onFormat = useCallback((next: CodeFormat) => update({ codeFormat: next }), [update]);
  const closeImport = useCallback(() => setImporting(false), []);

  return (
    <Panel
      title="Command"
      actions={
        <>
          <OptionToggle label="Command format" value={format} onChange={onFormat} options={CODE_FORMATS} />
          <CopyButton
            text={text}
            label={format === "cli" ? "Copy the command" : "Copy the config"}
            disabled={code === null}
          />
        </>
      }
    >
      <Collapsible open={importing} onOpenChange={setImporting} className="grid gap-2.5 p-3.5">
        <pre
          aria-label={format === "cli" ? "Command line" : "YAML config"}
          className="bg-muted/50 max-h-80 overflow-auto overscroll-contain rounded-md border px-3 py-2.5 font-mono text-xs leading-relaxed"
        >
          {code === null && (
            <span className="text-muted-foreground block">The command appears when the settings are valid.</span>
          )}
          {keyed.map(({ key, line, index }) => (
            <CodeLineView
              key={key}
              line={line}
              continued={format === "cli" && index < lines.length - 1}
              indent={format === "cli" && index > 0}
            />
          ))}
        </pre>
        <div className="flex flex-wrap items-center gap-2">
          <p className="text-muted-foreground flex-1 text-xs">
            {format === "cli"
              ? "Inputs, the settings that differ from the defaults, and the outputs the app adds to every run."
              : code === null
                ? ""
                : `Save as ${code.config_file} and run ${code.config_command}.`}
          </p>
          <CollapsibleTrigger render={<Button type="button" variant="outline" size="sm" />}>
            <FileUp aria-hidden />
            Load YAML
          </CollapsibleTrigger>
        </div>
        <CollapsibleContent>
          <YamlImport command={command} close={closeImport} />
        </CollapsibleContent>
      </Collapsible>
    </Panel>
  );
}

export function keyedLines(lines: readonly CodeLine[]): Array<{ key: string; line: CodeLine; index: number }> {
  const seen = new Map<string, number>();

  return lines.map((line, index) => {
    const base = `${line.kind}:${line.text}`;
    const count = seen.get(base) ?? 0;

    seen.set(base, count + 1);

    return { key: `${base}:${count}`, line, index };
  });
}

export function CodeLineView({ line, continued, indent }: { line: CodeLine; continued: boolean; indent: boolean }) {
  const prefix = indent && line.kind !== "comment" ? "  " : "";
  const suffix = continued && line.kind !== "comment" ? " \\" : "";

  return (
    <span
      className={cn(
        "block whitespace-pre",
        line.kind === "changed" && "text-primary font-bold",
        line.kind === "comment" && "text-muted-foreground",
      )}
    >
      {`${prefix}${line.text}${suffix}`}
    </span>
  );
}

function YamlImport({ command, close }: { command: AppCommand; close: () => void }) {
  const [text, setText] = useState("");
  const [messages, setMessages] = useState<string[]>([]);
  const [busy, setBusy] = useState(false);
  const loadConfig = useConfigLoader();
  const toasts = useToastManager();

  const apply = useCallback(async () => {
    setBusy(true);

    try {
      const result = await loadConfig(text, command, true);

      if (result.loaded) {
        const note = commandSwitchNote(command, result.command);

        if (note !== undefined) {
          toasts.add({ title: note });
        }

        close();
      } else {
        setMessages(result.messages);
      }
    } catch (error: unknown) {
      setMessages([errorMessage(error)]);
    } finally {
      setBusy(false);
    }
  }, [close, command, loadConfig, text, toasts]);

  const onApply = useCallback(() => void apply(), [apply]);
  const onText = useCallback((event: React.ChangeEvent<HTMLTextAreaElement>) => setText(event.target.value), []);

  return (
    <div className="grid gap-1.5">
      <Textarea
        rows={6}
        value={text}
        onChange={onText}
        aria-label="YAML config"
        placeholder={"clock_rate: 0.0008\nconfidence: true"}
        className="font-mono text-xs whitespace-pre"
      />
      {messages.length > 0 && (
        <ul className="text-destructive grid gap-0.5 text-xs">
          {messages.map((message) => (
            <li key={message}>{message}</li>
          ))}
        </ul>
      )}
      <div className="flex items-center gap-2">
        <Button type="button" variant="outline" size="sm" onClick={onApply} disabled={busy || text.trim() === ""}>
          {busy && <Spinner />}
          Apply settings
        </Button>
        <span className="text-muted-foreground text-xs">or drop a .yaml file anywhere on the page</span>
      </div>
    </div>
  );
}
