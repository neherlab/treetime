import type { LocalFiles, WorkspaceShell } from "@neherlab/app-contracts";
import { createContext, useContext } from "react";

export const LocalFilesContext = createContext<LocalFiles | null>(null);

export const WorkspaceShellContext = createContext<WorkspaceShell | null>(null);

export function useLocalFiles(): LocalFiles | null {
  return useContext(LocalFilesContext);
}

export function useWorkspaceShell(): WorkspaceShell | null {
  return useContext(WorkspaceShellContext);
}
