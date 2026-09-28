import type { LocalFiles } from "@neherlab/app-contracts";
import { RouterProvider } from "@tanstack/react-router";

import "./ui/fonts";
import { LocalFilesContext } from "./platform";
import { router } from "./router";
import { Toaster } from "./ui/toast";

export interface AppProps {
  localFiles?: LocalFiles | undefined;
}

export function App({ localFiles }: AppProps) {
  return (
    <LocalFilesContext.Provider value={localFiles ?? null}>
      <Toaster>
        <RouterProvider router={router} />
      </Toaster>
    </LocalFilesContext.Provider>
  );
}
