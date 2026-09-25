import type { LocalFiles } from "@neherlab/app-contracts";
import { RouterProvider } from "@tanstack/react-router";

import "./ui/fonts";
import { LocalFilesContext } from "./platform";
import { router } from "./router";
import { Toast } from "./ui";

export interface AppProps {
  localFiles?: LocalFiles | undefined;
}

export function App({ localFiles }: AppProps) {
  return (
    <LocalFilesContext.Provider value={localFiles ?? null}>
      <Toast.Provider>
        <RouterProvider router={router} />
        <Toast.Viewport />
      </Toast.Provider>
    </LocalFilesContext.Provider>
  );
}
