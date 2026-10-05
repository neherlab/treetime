import { RouterProvider } from "@tanstack/react-router";

import "./ui/fonts";
import type { Host } from "./host";
import { HostContext } from "./host-context";
import { router } from "./router";
import { Toaster } from "./ui/toast";

export interface AppProps {
  host: Host | null;
}

export function App({ host }: AppProps) {
  return (
    <HostContext.Provider value={host}>
      <Toaster>
        <RouterProvider router={router} />
      </Toaster>
    </HostContext.Provider>
  );
}
