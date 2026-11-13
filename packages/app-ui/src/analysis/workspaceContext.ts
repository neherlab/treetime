import type { AppCommand, CheckConfigResponse, InputFacts, SparseConfig } from "@neherlab/app-contracts";
import { createContext, use } from "react";

import type { RunBlock } from "./runAction";

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
  settled: SparseConfig;
  facts: InputFacts | undefined;
  check: CheckConfigResponse | undefined;
  checkPending: boolean;
  run: () => void;
  runBlock: RunBlock | undefined;
}
