import type { CodeLine, ConfigCode, UiCodeFormat } from "@neherlab/app-contracts";
import { Alert, AlertDescription } from "@neherlab/design/ui/alert";
import { Button } from "@neherlab/design/ui/button";
import { CopyButton } from "@neherlab/design/ui/copy-button";
import { Dialog, DialogContent, DialogDescription, DialogHeader, DialogTitle } from "@neherlab/design/ui/dialog";
import { Spinner } from "@neherlab/design/ui/spinner";
import { cn } from "cn";
import { useCallback, useMemo, type ReactNode } from "react";
import Download from "~icons/lucide/download";

import { OptionToggle } from "../components/OptionToggle";
import { saveText } from "../runs/saveRunOutput";
import { useDraftStore } from "../store/draft";

const CODE_FORMATS: ReadonlyArray<{ value: UiCodeFormat; label: string }> = [
  { value: "cli", label: "CLI" },
  { value: "yaml", label: "YAML" },
];

export type CommandLineCode = { state: "loading" } | { state: "invalid"; message: string } | ConfigCodeReady;

interface ConfigCodeReady {
  state: "ready";
  code: ConfigCode;
}

export function CommandLineDialog({
  open,
  onOpenChange,
  description,
  code,
}: {
  open: boolean;
  onOpenChange: (open: boolean) => void;
  description: ReactNode;
  code: CommandLineCode;
}) {
  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="max-h-[calc(100%-4rem)] grid-rows-[auto_minmax(0,1fr)] gap-4 sm:max-w-3xl">
        <DialogHeader>
          <DialogTitle>Command line</DialogTitle>
          <DialogDescription>{description}</DialogDescription>
        </DialogHeader>
        <div className="grid min-h-0 content-start gap-5 overflow-y-auto overscroll-contain">
          {code.state === "loading" && (
            <p className="text-muted-foreground flex items-center gap-2">
              <Spinner />
              Resolving the settings
            </p>
          )}
          {code.state === "invalid" && (
            <Alert variant="destructive">
              <AlertDescription>The command appears when the settings are valid: {code.message}</AlertDescription>
            </Alert>
          )}
          {code.state === "ready" && <CodeSections code={code.code} />}
        </div>
      </DialogContent>
    </Dialog>
  );
}

function CodeSections({ code }: { code: ConfigCode }) {
  const format = useDraftStore((state) => state.draft.code_format);
  const update = useDraftStore((state) => state.update);
  const onFormat = useCallback((next: UiCodeFormat) => update({ code_format: next }), [update]);

  const lines = format === "cli" ? code.command_line : code.yaml;
  const text = format === "cli" ? code.command_line_text : code.yaml_text;
  const keyed = useMemo(() => keyedLines(lines), [lines]);
  const citation = code.citation;

  const onSaveYaml = useCallback(
    () => saveText(code.yaml_text, code.config_file, "application/yaml"),
    [code.config_file, code.yaml_text],
  );

  return (
    <>
      <section aria-label="Command" className="grid gap-2">
        <div className="flex flex-wrap items-center gap-2">
          <OptionToggle label="Command format" value={format} onChange={onFormat} options={CODE_FORMATS} />
          <CopyButton text={text} label={format === "cli" ? "Copy the command" : "Copy the config"} />
          {format === "yaml" && (
            <Button type="button" variant="outline" size="sm" className="ml-auto" onClick={onSaveYaml}>
              <Download aria-hidden />
              Save {code.config_file}
            </Button>
          )}
        </div>
        <pre
          aria-label={format === "cli" ? "Command line" : "YAML config"}
          className="bg-muted/50 max-h-80 overflow-auto overscroll-contain rounded-md border px-3 py-2.5 font-mono text-xs leading-relaxed"
        >
          {keyed.map(({ key, line, index }) => (
            <CodeLineView
              key={key}
              line={line}
              continued={format === "cli" && index < lines.length - 1}
              indent={format === "cli" && index > 0}
            />
          ))}
        </pre>
        <p className="text-muted-foreground text-xs">
          {format === "cli" ? (
            "Inputs, the settings that differ from the defaults, and the outputs that the app adds to every run."
          ) : (
            <>
              Save as <code className="font-mono">{code.config_file}</code> and run{" "}
              <code className="font-mono">{code.config_command}</code>.
            </>
          )}
        </p>
      </section>
      <TextSection title="Methods" text={code.methods} copyLabel="Copy the methods paragraph">
        <p className="leading-relaxed">{code.methods}</p>
      </TextSection>
      <TextSection title="Citation" text={citation.text} copyLabel="Copy the citation">
        <p className="leading-relaxed">
          {citation.text}{" "}
          <a
            href={citation.url}
            target="_blank"
            rel="noopener noreferrer"
            className="text-primary font-bold underline-offset-4 hover:underline"
          >
            doi:{citation.doi}
          </a>
        </p>
      </TextSection>
    </>
  );
}

function keyedLines(lines: readonly CodeLine[]): Array<{ key: string; line: CodeLine; index: number }> {
  const seen = new Map<string, number>();

  return lines.map((line, index) => {
    const base = `${line.kind}:${line.text}`;
    const count = seen.get(base) ?? 0;

    seen.set(base, count + 1);

    return { key: `${base}:${count}`, line, index };
  });
}

function TextSection({
  title,
  text,
  copyLabel,
  children,
}: {
  title: string;
  text: string;
  copyLabel: string;
  children: ReactNode;
}) {
  return (
    <section aria-label={title} className="grid gap-1.5">
      <div className="flex items-center gap-2 border-b pb-1">
        <h3 className="font-heading text-sm font-bold">{title}</h3>
        <CopyButton text={text} label={copyLabel} />
      </div>
      {children}
    </section>
  );
}

function CodeLineView({ line, continued, indent }: { line: CodeLine; continued: boolean; indent: boolean }) {
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
