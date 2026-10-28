import { Toaster } from "@neherlab/design/ui/toast";
import "@neherlab/design/fonts";
import { RouterProvider } from "@tanstack/react-router";

import { type BrowserBackend, BrowserBackendContext } from "./browser-backend";
import type { Host } from "./host";
import { HostContext } from "./host-context";
import { LocaleProvider } from "./i18n/locale";
import { router } from "./router";

export interface AppProps {
  host: Host | null;
  browserBackend?: BrowserBackend;
}

export function App({ host, browserBackend }: AppProps) {
  return (
    <LocaleProvider locale="en">
      <HostContext.Provider value={host}>
        <BrowserBackendContext.Provider value={browserBackend ?? null}>
          <Toaster>
            <RouterProvider router={router} />
          </Toaster>
        </BrowserBackendContext.Provider>
      </HostContext.Provider>
    </LocaleProvider>
  );
}
