import type { RunRecord, RunResults, RunSummary } from "@neherlab/app-contracts";
import { runsAuspice, runsList, runsResults } from "@neherlab/app-contracts/client";
import { createContext, use, useMemo } from "react";

import { useApi } from "../api/hooks";
import type { RunProgress } from "../results/progress";
import type { TreeData } from "./TreeView";

export interface RunView {
  record: RunRecord;
  progress: RunProgress;
  failure: string | undefined;
}

export const RunContext = createContext<RunView | undefined>(undefined);

export function useRunView(): RunView {
  const view = use(RunContext);

  if (view === undefined) {
    throw new Error("A run view renders only inside the run page");
  }

  return view;
}

export function useRunSummary(id: string): RunSummary | undefined {
  const { data: list } = useApi((context) => runsList(context));

  return list?.runs.find((run) => run.id === id);
}

export type RunOutputs =
  | { state: "loading" }
  | { state: "failed"; message: string }
  | { state: "ready"; results: RunResults; tree: TreeData | undefined };

export function useRunOutputs(record: RunRecord): RunOutputs {
  const { data: results, error } = useApi((context) => runsResults({ ...context, path: { id: record.id } }), {
    staleTime: Infinity,
  });

  const hasTree = results?.tree !== undefined;

  const { data: document, error: treeError } = useApi(
    (context) => runsAuspice({ ...context, path: { id: record.id } }),
    { enabled: hasTree, staleTime: Infinity },
  );

  return useMemo((): RunOutputs => {
    const failure = error ?? treeError;

    if (failure !== null) {
      return { state: "failed", message: failure.message };
    }

    if (results === undefined || (hasTree && document === undefined)) {
      return { state: "loading" };
    }

    const tree =
      results.tree === undefined || document === undefined
        ? undefined
        : { document, tree: results.tree, run: { id: record.id, title: record.title, command: record.command } };

    return { state: "ready", results, tree };
  }, [document, error, hasTree, record.command, record.id, record.title, results, treeError]);
}
