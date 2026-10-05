import { createApiClient } from "@neherlab/app-contracts/client";
import {
  ApiProvider,
  apiPreferencesStorage,
  App,
  ErrorBoundary,
  PreferencesProvider,
  QueryProvider,
  ThemeProvider,
} from "@neherlab/app-ui";
import type { Host } from "@neherlab/app-ui/host";
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";

import { windowFetchConnection } from "../src/ipc-renderer";
import { createPortFetch } from "../src/port-fetch";

import "./index.css";

declare global {
  interface Window {
    treetimeHost: Host;
  }
}

const host = window.treetimeHost;

const connection = windowFetchConnection(window, host);

const client = createApiClient({ baseUrl: globalThis.location.origin, fetch: createPortFetch(connection) });

const preferences = apiPreferencesStorage(client);

const root = document.getElementById("root");

if (root) {
  createRoot(root).render(
    <StrictMode>
      <ThemeProvider>
        <ErrorBoundary>
          <QueryProvider>
            <ApiProvider client={client}>
              <PreferencesProvider storage={preferences}>
                <App host={host} />
              </PreferencesProvider>
            </ApiProvider>
          </QueryProvider>
        </ErrorBoundary>
      </ThemeProvider>
    </StrictMode>,
  );
}
