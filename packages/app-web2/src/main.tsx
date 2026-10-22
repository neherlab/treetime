import { createApiClient } from "@neherlab/app-contracts/client";
import {
  ApiProvider,
  App,
  browserPreferencesStorage,
  ErrorBoundary,
  PreferencesProvider,
  QueryProvider,
  reloadOnChunkError,
  TAB_BACKEND,
  ThemeProvider,
} from "@neherlab/app-ui";
import { createPortFetch } from "@neherlab/app-ui/port-fetch";
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";

import { browserWorkerHost } from "./backend/browser-host";
import { WorkerConnection } from "./backend/connection";
import { parseJobs } from "./backend/threads";

import "./index.css";

reloadOnChunkError(globalThis.window);

const { jobs, warning } = parseJobs(globalThis.location.search);

if (warning !== undefined) {
  console.warn(`[TreeTime] ${warning}`);
}

const connection = new WorkerConnection(browserWorkerHost(jobs));

const client = createApiClient({
  baseUrl: globalThis.location.origin,
  fetch: createPortFetch(connection, { transferBodies: true }),
});

const preferences = browserPreferencesStorage(globalThis.localStorage, globalThis.sessionStorage);

const root = document.getElementById("root");

if (root) {
  createRoot(root).render(
    <StrictMode>
      <ThemeProvider>
        <ErrorBoundary>
          <QueryProvider>
            <ApiProvider client={client}>
              <PreferencesProvider storage={preferences}>
                <App host={null} browserBackend={TAB_BACKEND} />
              </PreferencesProvider>
            </ApiProvider>
          </QueryProvider>
        </ErrorBoundary>
      </ThemeProvider>
    </StrictMode>,
  );
}
