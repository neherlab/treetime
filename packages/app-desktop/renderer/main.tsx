import { createApiClient } from "@neherlab/app-contracts/client";
import { ApiProvider, App, ErrorBoundary, QueryProvider, ThemeProvider } from "@neherlab/app-ui";
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";

import {
  createDesktopSaveActions,
  createLocalFiles,
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

const DESKTOP_ORIGIN = "http://treetime.desktop";

const connection = windowFetchConnection(window, window.treetimeShell);

const client = createApiClient({ baseUrl: DESKTOP_ORIGIN, fetch: createPortFetch(connection) });

const save = createDesktopSaveActions(window.treetimeShell);

const localFiles = createLocalFiles(window.treetimeShell);

const root = document.getElementById("root");

if (root) {
  createRoot(root).render(
    <StrictMode>
      <ThemeProvider>
        <ErrorBoundary>
          <QueryProvider>
            <ApiProvider client={client} save={save}>
              <App localFiles={localFiles} />
            </ApiProvider>
          </QueryProvider>
        </ErrorBoundary>
      </ThemeProvider>
    </StrictMode>,
  );
}
