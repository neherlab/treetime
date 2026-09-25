import type { LogEvent, RunRecordResult } from "@neherlab/app-contracts";
import { useQueryClient } from "@tanstack/react-query";
import { Link, useNavigate } from "@tanstack/react-router";
import { DateTime } from "luxon";
import { useCallback, useEffect, useState } from "react";

import { defaultText } from "../analysis/SettingField";
import { useBridge } from "../BridgeContext";
import { formatDuration } from "../format";
import { RUNS_KEY, useRunRecord } from "../queries";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { COMMAND_INFO } from "../settings/commands";
import { settingValue } from "../settings/config";
import { zJsonObject } from "../settings/json";
import { settingLabel } from "../settings/labels";
import { rerunDraft } from "../settings/rerun";
import { StatusIcon, statusLabel } from "../shell/StatusIcon";
import { useDraftStore } from "../store/draft";
import { Button, Toast, cn } from "../ui";

export type RunTab = "results" | "settings" | "log";

const TABS: ReadonlyArray<{
  tab: RunTab;
  label: string;
  to: "/runs/$id/results" | "/runs/$id/settings" | "/runs/$id/log";
}> = [
  { tab: "results", label: "Results", to: "/runs/$id/results" },
  { tab: "settings", label: "Settings", to: "/runs/$id/settings" },
  { tab: "log", label: "Log", to: "/runs/$id/log" },
];

export function RunPage({ id, tab }: { id: string; tab: RunTab }) {
  const { data: record, error } = useRunRecord(id);

  if (error !== null) {
    return <p className="text-signal-danger p-10 text-center">The run cannot be loaded: {error.message}</p>;
  }

  if (record === undefined) {
    return <p className="text-ink-muted p-10 text-center">Loading the run...</p>;
  }

  return (
    <div className="mx-auto max-w-[92.5rem] px-5 pt-4 pb-16">
      <RunHeader record={record} />
      <nav aria-label="Run views" className="border-line mb-4 flex gap-0.5 border-b">
        {TABS.map((entry) => (
          <Link
            key={entry.tab}
            to={entry.to}
            params={{ id }}
            aria-current={entry.tab === tab ? "page" : undefined}
            className="text-ink-muted aria-[current=page]:border-accent aria-[current=page]:text-ink -mb-px border-b-2 border-transparent px-3 py-2 font-bold"
          >
            {entry.label}
          </Link>
        ))}
      </nav>
      {tab === "results" && <RunSummaryView record={record} />}
      {tab === "settings" && <RunSettingsView record={record} />}
      {tab === "log" && <RunLogView id={id} />}
    </div>
  );
}

function RunHeader({ record }: { record: RunRecordResult }) {
  const bridge = useBridge();
  const navigate = useNavigate();
  const queryClient = useQueryClient();
  const toasts = Toast.useToastManager();
  const created = DateTime.fromISO(record.created_at);

  const editAndRunAgain = useCallback(() => {
    const draft = rerunDraft(record);

    useDraftStore.getState().load({
      command: record.command,
      config: draft.config,
      sources: Object.fromEntries(
        Object.entries(draft.inputLabels).map(([key, label]) => [key, { label, origin: "run" as const, size: null }]),
      ),
      title: draft.title,
      fromRunId: record.id,
      uploadRunId: null,
    });
    void navigate({ to: "/new" });
  }, [navigate, record]);

  const cancel = useCallback(async () => {
    try {
      await bridge.cancelRun(record.id);
      await queryClient.invalidateQueries({ queryKey: RUNS_KEY });
    } catch (error: unknown) {
      toasts.add({ title: "The run cannot be cancelled", description: error instanceof Error ? error.message : "" });
    }
  }, [bridge, queryClient, record.id, toasts]);

  const onCancel = useCallback(() => void cancel(), [cancel]);

  const togglePin = useCallback(async () => {
    try {
      await bridge.updateRun(record.id, { pinned: !record.pinned });
      await queryClient.invalidateQueries({ queryKey: RUNS_KEY });
    } catch (error: unknown) {
      toasts.add({ title: "The run cannot be pinned", description: error instanceof Error ? error.message : "" });
    }
  }, [bridge, queryClient, record.id, record.pinned, toasts]);

  const onTogglePin = useCallback(() => void togglePin(), [togglePin]);

  return (
    <div className="mb-3.5 flex flex-wrap items-start gap-4">
      <div className="min-w-0">
        <h1 className="text-2xl leading-tight font-bold">{record.title}</h1>
        <div className="text-ink-muted mt-1 flex flex-wrap items-center gap-x-3 gap-y-1.5">
          <span className="inline-flex items-center gap-1.5 font-bold">
            <StatusIcon status={record.status} />
            {statusLabel(record.status)}
          </span>
          <span>{COMMAND_INFO[record.command].label}</span>
          <span>{created.isValid ? created.toFormat("d LLL yyyy, HH:mm") : record.created_at}</span>
          {record.duration_seconds !== null && record.duration_seconds !== undefined && (
            <span>{formatDuration(record.duration_seconds)}</span>
          )}
        </div>
      </div>
      <div className="ml-auto flex gap-1.5">
        {record.status === "running" && (
          <Button type="button" variant="outline" size="sm" onClick={onCancel}>
            Cancel
          </Button>
        )}
        <Button type="button" variant="ghost" size="sm" onClick={onTogglePin}>
          {record.pinned ? "Unpin" : "Pin"}
        </Button>
        <Button type="button" variant="outline" size="sm" onClick={editAndRunAgain}>
          Edit and run again
        </Button>
      </div>
    </div>
  );
}

function RunSummaryView({ record }: { record: RunRecordResult }) {
  const error = record.error;

  return (
    <div className="grid gap-3.5">
      {error !== null && error !== undefined && (
        <div className="border-signal-danger bg-signal-danger-subtle rounded-lg border px-4 py-3.5">
          <h3 className="mb-1.5 text-[0.9375rem] font-bold">The run failed</h3>
          <p>{error.message}</p>
          {error.causes.length > 0 && (
            <ul className="text-ink-muted mt-1 list-disc pl-5">
              {error.causes.map((cause) => (
                <li key={cause}>{cause}</li>
              ))}
            </ul>
          )}
        </div>
      )}
      <div className="border-line bg-surface-1 rounded-lg border">
        <h3 className="border-line border-b px-3.5 py-2.5 font-bold">Output files</h3>
        {record.output_files.length === 0 ? (
          <p className="text-ink-muted px-3.5 py-3">
            {record.status === "running" ? "The run is writing its outputs." : "The run wrote no files."}
          </p>
        ) : (
          <ul className="px-3.5 py-2">
            {record.output_files.map((file) => (
              <li key={file.path} className="border-line flex justify-between gap-3 border-b py-1.5 last:border-b-0">
                <code className="font-mono text-xs">{file.path}</code>
                <span className="text-ink-faint text-xs">{file.kind}</span>
              </li>
            ))}
          </ul>
        )}
      </div>
    </div>
  );
}

function RunSettingsView({ record }: { record: RunRecordResult }) {
  const specs = COMMAND_SETTINGS[record.command].specs;
  const config = zJsonObject.parse(record.config);
  const changed = new Set(record.changed_settings);
  const rows = specs.filter((spec) => changed.has(spec.key));

  return (
    <div className="border-line bg-surface-1 rounded-lg border">
      <h3 className="border-line border-b px-3.5 py-2.5 font-bold">Settings that differ from the defaults</h3>
      {rows.length === 0 ? (
        <p className="text-ink-muted px-3.5 py-3">Every setting has its default value.</p>
      ) : (
        <table className="w-full border-collapse">
          <thead>
            <tr className="text-ink-faint text-left text-xs">
              <th className="px-3.5 py-1.5">Setting</th>
              <th className="px-3.5 py-1.5">Value</th>
              <th className="px-3.5 py-1.5">Default</th>
            </tr>
          </thead>
          <tbody>
            {rows.map((spec) => (
              <tr key={spec.key} className="border-line border-t">
                <td className="px-3.5 py-1.5">
                  {settingLabel(spec.key)} <code className="text-ink-faint font-mono text-xs">{spec.flag}</code>
                </td>
                <td className="px-3.5 py-1.5 font-mono text-xs">{defaultText(settingValue(config, spec))}</td>
                <td className="text-ink-faint px-3.5 py-1.5 font-mono text-xs">{defaultText(spec.defaultValue)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      )}
    </div>
  );
}

interface LogLine {
  seq: number;
  event: LogEvent;
}

function RunLogView({ id }: { id: string }) {
  const bridge = useBridge();
  const [lines, setLines] = useState<LogLine[]>([]);
  const [failure, setFailure] = useState<string | null>(null);

  useEffect(() => {
    const controller = new AbortController();
    const received: LogLine[] = [];

    bridge
      .followRun(id, {
        signal: controller.signal,
        onEvent: (event) => {
          if (event.type === "log") {
            received.push({ seq: event.seq, event: event.data });
            setLines([...received]);
          }
        },
      })
      .catch((error: Error) => {
        if (!controller.signal.aborted) {
          setFailure(error.message);
        }
      });

    return () => controller.abort();
  }, [bridge, id]);

  return (
    <div className="border-line bg-surface-2 max-h-[70vh] overflow-auto rounded-md border px-2.5 py-2 font-mono text-[0.71875rem] leading-relaxed whitespace-pre-wrap">
      {failure !== null && <p className="text-signal-danger">{failure}</p>}
      {lines.length === 0 && failure === null && <p className="text-ink-faint">No log lines yet.</p>}
      {lines.map((line) => (
        <div
          key={line.seq}
          className={cn(
            line.event.level === "warn" && "text-signal-warn",
            line.event.level === "error" && "text-signal-danger font-bold",
          )}
        >
          {line.event.message}
        </div>
      ))}
    </div>
  );
}
