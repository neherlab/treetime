import { ApiProvider, App, ErrorBoundary, QueryProvider, reloadOnChunkError, ThemeProvider } from "@neherlab/app-ui";
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";

import { createWebApiClient } from "./api-client";
import { createWebSaveActions } from "./save-web";

import "./index.css";

reloadOnChunkError(globalThis.window);

const client = createWebApiClient();

const save = createWebSaveActions(client);

const root = document.getElementById("root");

if (root) {
  createRoot(root).render(
    <StrictMode>
      <ThemeProvider>
        <ErrorBoundary>
          <QueryProvider>
            <ApiProvider client={client} save={save}>
              <App />
            </ApiProvider>
          </QueryProvider>
        </ErrorBoundary>
      </ThemeProvider>
    </StrictMode>,
  );
}
