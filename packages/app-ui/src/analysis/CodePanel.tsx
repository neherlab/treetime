import { errorMessage } from "@neherlab/app-contracts";
import type { AppCommand, CodeLine, ConfigCode } from "@neherlab/app-contracts";
import { Copy } from "lucide-react";
import { useCallback, useMemo, useState } from "react";

import { commandSwitchNote } from "../settings/commands";
import { useDraftStore } from "../store/draft";
import type { CodeFormat } from "../store/draftSchema";
import { Button, Segmented, Toast, cn } from "../ui";
import { useConfigLoader } from "./useConfigLoader";

const CODE_FORMATS: ReadonlyArray<{ value: CodeFormat; label: string }> = [
  { value: "cli", label: "CLI" },
  { value: "yaml", label: "YAML" },
];

export function CodePanel({ command, code }: { command: AppCommand; code: ConfigCode | null }) {
  const format = useDraftStore((state) => state.codeFormat);
  const update = useDraftStore((state) => state.update);
  const toasts = Toast.useToastManager();
  const [importing, setImporting] = useState(false);

  const lines = useMemo(() => (code === null ? [] : format === "cli" ? code.command_line : code.yaml), [code, format]);
  const text = code === null ? "" : format === "cli" ? code.command_line_text : code.yaml_text;
  const keyed = useMemo(() => keyedLines(lines), [lines]);

  const onFormat = useCallback((next: CodeFormat) => update({ codeFormat: next }), [update]);
  const toggleImport = useCallback(() => setImporting((shown) => !shown), []);
  const closeImport = useCallback(() => setImporting(false), []);

  const copy = useCallback(() => {
    navigator.clipboard.writeText(text).then(
      () =>
        toasts.add({ title: format === "cli" ? "Command copied to the clipboard" : "Config copied to the clipboard" }),
      () => toasts.add({ title: "The browser blocked clipboard access" }),
    );
  }, [format, text, toasts]);

  return (
    <div className="border-line bg-surface-1 rounded-lg border">
      <div className="border-line flex items-center gap-2.5 border-b px-3.5 py-2.5">
        <h3 className="font-bold">Command</h3>
        <div className="ml-auto flex items-center gap-2">
          <Segmented label="Command format" value={format} onChange={onFormat} options={CODE_FORMATS} />
          <Button type="button" variant="ghost" size="icon" aria-label="Copy" onClick={copy} disabled={code === null}>
            <Copy size={14} aria-hidden />
          </Button>
        </div>
      </div>
      <div className="grid gap-2.5 px-3.5 py-3">
        <pre
          aria-label={format === "cli" ? "Command line" : "YAML config"}
          className="border-line bg-surface-2 m-0 max-h-80 overflow-auto rounded-md border px-3 py-2.5 font-mono text-xs leading-relaxed"
        >
          {code === null && (
            <span className="text-ink-faint block">The command appears when the settings are valid.</span>
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
        <p className="text-ink-faint text-xs">
          {format === "cli"
            ? "Inputs, the settings that differ from the defaults, and the outputs the app adds to every run. "
            : code === null
              ? ""
              : `Save as ${code.config_file} and run ${code.config_command}. `}
          <button type="button" onClick={toggleImport} className="text-accent font-bold">
            Load YAML
          </button>
        </p>
        {importing && <YamlImport command={command} close={closeImport} />}
      </div>
    </div>
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
        line.kind === "changed" && "text-accent font-bold",
        line.kind === "comment" && "text-ink-faint",
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
  const toasts = Toast.useToastManager();

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
      <textarea
        rows={6}
        value={text}
        onChange={onText}
        aria-label="YAML config"
        placeholder={"clock_rate: 0.0008\nconfidence: true"}
        className="border-line bg-surface-2 rounded-md border px-3 py-2 font-mono text-xs whitespace-pre"
      />
      {messages.length > 0 && (
        <ul className="text-signal-danger grid gap-0.5 text-xs">
          {messages.map((message) => (
            <li key={message}>{message}</li>
          ))}
        </ul>
      )}
      <div className="flex items-center gap-2">
        <Button type="button" variant="outline" size="sm" onClick={onApply} disabled={busy || text.trim() === ""}>
          Apply settings
        </Button>
        <span className="text-ink-faint text-xs">or drop a .yaml file anywhere on the page</span>
      </div>
    </div>
  );
}
