import { createApiClient } from "@neherlab/app-contracts/client";
import {
  ApiProvider,
  App,
  browserPreferencesStorage,
  ErrorBoundary,
  PreferencesProvider,
  QueryProvider,
  reloadOnChunkError,
  ThemeProvider,
} from "@neherlab/app-ui";
import { createPortFetch } from "@neherlab/app-ui/port-fetch";
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";

import { browserWorkerHost } from "./backend/browser-host";
import { WorkerConnection } from "./backend/connection";

import "./index.css";

reloadOnChunkError(globalThis.window);

const connection = new WorkerConnection(browserWorkerHost);

const client = createApiClient({ baseUrl: globalThis.location.origin, fetch: createPortFetch(connection) });

const preferences = browserPreferencesStorage(globalThis.localStorage);

const root = document.getElementById("root");

if (root) {
  createRoot(root).render(
    <StrictMode>
      <ThemeProvider>
        <ErrorBoundary>
          <QueryProvider>
            <ApiProvider client={client}>
              <PreferencesProvider storage={preferences}>
                <App host={null} />
              </PreferencesProvider>
            </ApiProvider>
          </QueryProvider>
        </ErrorBoundary>
      </ThemeProvider>
    </StrictMode>,
  );
}
