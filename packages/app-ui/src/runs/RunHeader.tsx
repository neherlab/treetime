import type { RunRecord } from "@neherlab/app-contracts";
import { runsList } from "@neherlab/app-contracts/client";
import { Button } from "@neherlab/design/ui/button";
import { NativeSelect, NativeSelectOption } from "@neherlab/design/ui/native-select";
import { useNavigate } from "@tanstack/react-router";
import { DateTime } from "luxon";
import { useCallback, useMemo } from "react";
import RotateCcw from "~icons/lucide/rotate-ccw";

import { useApi } from "../api/hooks";
import { formatDuration } from "../format";
import { commandSettings } from "../settings/catalog";
import { baseName } from "../settings/inputs";
import { RunActionsMenu } from "../shell/runActions";
import { StatusIcon, statusLabel } from "../shell/StatusIcon";
import { RunTitle } from "./RunTitle";
import { useRerun } from "./useRerun";

const TREE_SETTING = "tree";

export function RunHeader({ record }: { record: RunRecord }) {
  const created = DateTime.fromISO(record.created_at);
  const dataset = runDataset(record);

  return (
    <header className="flex flex-wrap items-start gap-x-4 gap-y-2">
      <div className="grid min-w-0 gap-1">
        <RunTitle record={record} />
        <p className="text-muted-foreground flex flex-wrap items-center gap-x-3 gap-y-1 text-sm">
          <span className="text-foreground inline-flex items-center gap-1.5 font-bold">
            <StatusIcon status={record.status} />
            {statusLabel(record.status)}
          </span>
          <span>{commandSettings(record.command).title}</span>
          {dataset !== undefined && <span className="font-mono text-xs">{dataset}</span>}
          <time dateTime={record.created_at}>
            {created.isValid ? created.toFormat("d LLL yyyy, HH:mm") : record.created_at}
          </time>
          {record.duration_seconds !== undefined && <span>{formatDuration(record.duration_seconds)}</span>}
        </p>
      </div>
      <RunHeaderActions record={record} />
    </header>
  );
}

function RunHeaderActions({ record }: { record: RunRecord }) {
  const navigate = useNavigate();
  const rerun = useRerun(record);
  const { data: list } = useApi((context) => runsList(context));

  const comparable = useMemo(
    () => (list?.runs ?? []).filter((run) => run.id !== record.id && run.status === "ok"),
    [list, record.id],
  );

  const onCompare = useCallback(
    (event: React.ChangeEvent<HTMLSelectElement>) => {
      const other = event.target.value;

      if (other !== "") {
        void navigate({ to: "/compare/$a/$b", params: { a: record.id, b: other } });
      }
    },
    [navigate, record.id],
  );

  return (
    <div data-slot="run-actions" className="ml-auto flex flex-wrap items-center gap-1.5">
      {comparable.length > 0 && (
        <NativeSelect size="sm" aria-label="Compare with another run" value="" onChange={onCompare} className="w-44">
          <NativeSelectOption value="">Compare with...</NativeSelectOption>
          {comparable.map((run) => (
            <NativeSelectOption key={run.id} value={run.id}>
              {run.title}
            </NativeSelectOption>
          ))}
        </NativeSelect>
      )}
      <Button type="button" variant="outline" size="sm" onClick={rerun}>
        <RotateCcw aria-hidden />
        Edit and run again
      </Button>
      <RunActionsMenu run={record} />
    </div>
  );
}

function runDataset(record: RunRecord): string | undefined {
  if (record.example !== undefined) {
    return record.example;
  }

  const tree = record.inputs.find((input) => input.setting === TREE_SETTING);

  return tree === undefined ? undefined : baseName(tree.path);
}
