import type { RunRecord } from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { DateTime } from "luxon";
import GitBranchPlus from "~icons/lucide/git-branch-plus";

import { formatDuration } from "../format";
import { commandSettings } from "../settings/catalog";
import { baseName } from "../settings/inputs";
import { hasEnded, RunActionsMenu } from "../shell/runActions";
import { StatusIcon, statusLabel } from "../shell/StatusIcon";
import { DownloadMenu } from "./DownloadMenu";
import { RunCommandLine } from "./RunCommandLine";
import { useRunSummary } from "./runContext";
import { RunTitle } from "./RunTitle";
import { useCreateVariant } from "./useCreateVariant";

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
  const createVariant = useCreateVariant(record);
  const summary = useRunSummary(record.id);

  return (
    <div data-slot="run-actions" className="ml-auto flex flex-wrap items-center gap-1.5">
      {hasEnded(record.status) && <DownloadMenu record={record} packageName={summary?.package_name} />}
      <RunCommandLine record={record} />
      <Button type="button" variant="outline" size="sm" onClick={createVariant}>
        <GitBranchPlus aria-hidden />
        Create variant
      </Button>
      {summary !== undefined && <RunActionsMenu run={summary} />}
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
