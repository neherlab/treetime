import { errorMessage, type JobId } from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { useCallback, useMemo } from "react";
import Play from "~icons/lucide/play";
import Undo from "~icons/lucide/undo-2";

import { useDraftStore } from "../store/draft";
import { ownDraft, runOptionsDraft } from "../store/draftRules";
import { RUN_BLOCK_TEXT } from "./runAction";
import { useShownRunSettings } from "./shownRun";
import { staleRun, staleText } from "./staleRun";
import { useWorkspace } from "./workspaceContext";

export function StaleBar({ shownRunId }: { shownRunId: JobId }) {
  const { check, run, runBlock } = useWorkspace();
  const edited = useDraftStore((state) => state.viewing === undefined);
  const { record, settings } = useShownRunSettings(shownRunId);

  const stale = useMemo(
    () => staleRun({ shownRunId, runCommand: record.data?.command, edited, check }),
    [check, edited, record.data?.command, shownRunId],
  );

  const revertible = record.data?.prepared === true && settings.data !== undefined;

  const revert = useCallback(() => {
    const store = useDraftStore.getState();

    const options =
      record.data === undefined || settings.data === undefined
        ? undefined
        : runOptionsDraft(record.data, settings.data, ownDraft(store));

    if (options !== undefined) {
      store.revertToRun(shownRunId, options);
    }
  }, [record.data, settings.data, shownRunId]);

  const loadError = record.error ?? settings.error;

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
          disabled={!revertible}
          title={
            loadError === null ? undefined : `The settings of the run cannot be loaded: ${errorMessage(loadError)}`
          }
        >
          <Undo aria-hidden />
          Revert to this run&apos;s settings
        </Button>
        <Button
          type="button"
          size="sm"
          onClick={run}
          disabled={runBlock !== undefined}
          title={runBlock === undefined ? undefined : RUN_BLOCK_TEXT[runBlock]}
        >
          <Play aria-hidden />
          Run
        </Button>
      </div>
      {loadError !== null && (
        <p role="alert" className="text-destructive basis-full text-xs">
          The settings of the run cannot be loaded: {errorMessage(loadError)}
        </p>
      )}
    </div>
  );
}
