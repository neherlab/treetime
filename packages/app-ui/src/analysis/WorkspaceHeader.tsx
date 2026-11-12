import type { CheckConfigResponse, JobId } from "@neherlab/app-contracts";
import { runsList } from "@neherlab/app-contracts/client";
import { Button } from "@neherlab/design/ui/button";
import { useCallback } from "react";
import RotateCcw from "~icons/lucide/rotate-ccw";

import { useApi } from "../api/hooks";
import { useDraftStore } from "../store/draft";
import { variantSummary } from "./variantSummary";

export function WorkspaceTitle({ check }: { check: CheckConfigResponse | undefined }) {
  return check?.dataset ?? "New analysis";
}

export function WorkspaceDescription({
  shownRunId,
  check,
}: {
  shownRunId: JobId | undefined;
  check: CheckConfigResponse | undefined;
}) {
  const viewing = useDraftStore((state) => state.viewing?.runId);
  const fromRunId = useDraftStore((state) => state.draft.from_run_id);
  const { data: runList } = useApi((context) => runsList(context));
  const findRun = (id: JobId | undefined) => runList?.runs.find((run) => run.id === id);

  const text =
    viewing === undefined
      ? variantSummary({ fromRunId, shownRunId, viewing: false, source: findRun(fromRunId), check })
      : `Options of "${findRun(viewing)?.title ?? viewing}". An edit starts a variant of this run.`;

  return text === undefined ? null : <p className="text-muted-foreground text-xs">{text}</p>;
}

export function ResetButton() {
  const command = useDraftStore((state) => state.draft.command);
  const discard = useDraftStore((state) => state.discard);
  const reset = useCallback(() => discard("reset", command), [command, discard]);

  return (
    <Button type="button" variant="ghost" size="sm" onClick={reset}>
      <RotateCcw aria-hidden />
      Reset form
    </Button>
  );
}
