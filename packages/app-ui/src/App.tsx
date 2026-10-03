import type { LocalFiles, WorkspaceShell } from "@neherlab/app-contracts";
import { RouterProvider } from "@tanstack/react-router";

import "./ui/fonts";
import { LocalFilesContext, WorkspaceShellContext } from "./platform";
import { router } from "./router";
import { Toaster } from "./ui/toast";

export interface AppProps {
  localFiles?: LocalFiles | undefined;
  workspaceShell?: WorkspaceShell | undefined;
}

export function App({ localFiles, workspaceShell }: AppProps) {
  return (
    <LocalFilesContext.Provider value={localFiles ?? null}>
      <WorkspaceShellContext.Provider value={workspaceShell ?? null}>
        <Toaster>
          <RouterProvider router={router} />
        </Toaster>
      </WorkspaceShellContext.Provider>
    </LocalFilesContext.Provider>
  );
}
