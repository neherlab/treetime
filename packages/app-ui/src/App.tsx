import { Toaster } from "@neherlab/design/ui/toast";
import "@neherlab/design/fonts";
import { RouterProvider } from "@tanstack/react-router";

import type { Host } from "./host";
import { HostContext } from "./host-context";
import { LocaleProvider } from "./i18n/locale";
import { router } from "./router";

export interface AppProps {
  host: Host | null;
}

export function App({ host }: AppProps) {
  return (
    <LocaleProvider locale="en">
      <HostContext.Provider value={host}>
        <Toaster>
          <RouterProvider router={router} />
        </Toaster>
      </HostContext.Provider>
    </LocaleProvider>
  );
}
