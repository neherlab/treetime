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
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";

import {
  createDesktopSaveActions,
  createLocalFiles,
  createWorkspaceShell,
  windowFetchConnection,
  type DesktopShell,
} from "../src/desktop-shell";
import { createPortFetch } from "../src/port-fetch";

import "./index.css";

declare global {
  interface Window {
    treetimeShell: DesktopShell;
  }
}

const connection = windowFetchConnection(window, window.treetimeShell);

const client = createApiClient({ baseUrl: globalThis.location.origin, fetch: createPortFetch(connection) });

const save = createDesktopSaveActions(window.treetimeShell);

const localFiles = createLocalFiles(window.treetimeShell);

const workspaceShell = createWorkspaceShell(window.treetimeShell);

const preferences = apiPreferencesStorage(client);

const root = document.getElementById("root");

if (root) {
  createRoot(root).render(
    <StrictMode>
      <ThemeProvider>
        <ErrorBoundary>
          <QueryProvider>
            <ApiProvider client={client} save={save}>
              <PreferencesProvider storage={preferences}>
                <App localFiles={localFiles} workspaceShell={workspaceShell} />
              </PreferencesProvider>
            </ApiProvider>
          </QueryProvider>
        </ErrorBoundary>
      </ThemeProvider>
    </StrictMode>,
  );
}
