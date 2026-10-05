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
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";

import { createWebApiClient } from "./api-client";

import "./index.css";

reloadOnChunkError(globalThis.window);

const client = createWebApiClient();

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
