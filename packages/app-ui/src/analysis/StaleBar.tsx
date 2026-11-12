import type { JobId } from "@neherlab/app-contracts";
import { runsGet, runsSettings } from "@neherlab/app-contracts/client";
import { Button } from "@neherlab/design/ui/button";
import { useCallback, useMemo } from "react";
import Play from "~icons/lucide/play";
import Undo from "~icons/lucide/undo-2";

import { useApi } from "../api/hooks";
import { commandSettings } from "../settings/catalog";
import { useDraftStore } from "../store/draft";
import { ownDraft, runOptionsDraft } from "../store/draftRules";
import { staleRun, staleText } from "./staleRun";
import { useWorkspace } from "./workspaceContext";

export function StaleBar({ shownRunId }: { shownRunId: JobId }) {
  const { command, check, run } = useWorkspace();
  const edited = useDraftStore((state) => state.viewing === undefined);

  const { data: record } = useApi((context) => runsGet({ ...context, path: { id: shownRunId } }));

  const { data: settings } = useApi((context) => runsSettings({ ...context, path: { id: shownRunId } }), {
    staleTime: Infinity,
  });

  const stale = useMemo(
    () =>
      staleRun({
        shownRunId,
        runCommand: record?.command,
        edited,
        check,
        specs: commandSettings(command).settings,
      }),
    [check, command, edited, record?.command, shownRunId],
  );

  const revert = useCallback(() => {
    if (record === undefined || settings === undefined) {
      return;
    }

    const store = useDraftStore.getState();

    store.revertToRun(shownRunId, runOptionsDraft(record, settings, ownDraft(store)));
  }, [record, settings, shownRunId]);

  const blocked = check?.status === "invalid" || (check?.checks.some((finding) => finding.level === "block") ?? true);

  if (!stale.stale) {
    return null;
  }

  return (
    <div className="bg-muted/60 sticky top-0 z-10 flex flex-wrap items-center gap-x-3 gap-y-1.5 border-b px-4 py-2 text-sm">
      <output className="min-w-0 flex-1">{staleText(stale.changes)}</output>
      <div className="flex flex-wrap items-center gap-1.5">
        <Button
          type="button"
          variant="ghost"
          size="sm"
          onClick={revert}
          disabled={record === undefined || settings === undefined}
        >
          <Undo aria-hidden />
          Revert to this run&apos;s settings
        </Button>
        <Button type="button" size="sm" onClick={run} disabled={blocked}>
          <Play aria-hidden />
          Run
        </Button>
      </div>
    </div>
  );
}
