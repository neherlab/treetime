import type { AppCommand, CheckConfigResponse, InputFacts, JobId, SparseConfig } from "@neherlab/app-contracts";
import { createContext, use } from "react";

export const WorkspaceContext = createContext<WorkspaceView | undefined>(undefined);

export function useWorkspace(): WorkspaceView {
  const view = use(WorkspaceContext);

  if (view === undefined) {
    throw new Error("useWorkspace must be used inside the dataset workspace");
  }

  return view;
}

export interface WorkspaceView {
  command: AppCommand;
  shownRunId: JobId | undefined;
  readOnly: boolean;
  settled: SparseConfig;
  facts: InputFacts | undefined;
  check: CheckConfigResponse | undefined;
  run: () => void;
}
