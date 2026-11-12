import type { AppCommand, RunSummary } from "@neherlab/app-contracts";
import { runsList } from "@neherlab/app-contracts/client";
import { useEffect, useMemo } from "react";

import { useApi } from "../api/hooks";
import { useHiddenRunIds, useRunDeletionStore } from "../store/runDeletion";
import { listedRuns } from "./runList";

export function useListedRuns(filter = "", command: AppCommand | null = null): ListedRuns {
  const { data, error } = useApi((context) => runsList(context));
  const hidden = useHiddenRunIds();
  const forget = useRunDeletionStore((state) => state.forget);

  useEffect(() => {
    if (data !== undefined) {
      forget(new Set(data.runs.map((run) => run.id)));
    }
  }, [data, forget]);

  const all = useMemo(() => listedRuns(data?.runs ?? [], "", null, hidden), [data, hidden]);
  const runs = useMemo(() => listedRuns(all, filter, command), [all, command, filter]);

  return { all, runs, loaded: data !== undefined, error };
}

interface ListedRuns {
  all: RunSummary[];
  runs: RunSummary[];
  loaded: boolean;
  error: Error | null;
}
