import { ApiProvider, App, BridgeProvider, ErrorBoundary, QueryProvider, ThemeProvider } from "@neherlab/app-ui";
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";

import { createWebApiClient } from "./api-client";
import { createWebBridge } from "./bridge-web";
import { createWebSaveActions } from "./save-web";

import "./index.css";

const client = createWebApiClient();

const bridge = createWebBridge({ client });

const save = createWebSaveActions(client);

const root = document.getElementById("root");

if (root) {
  createRoot(root).render(
    <StrictMode>
      <ThemeProvider>
        <ErrorBoundary>
          <BridgeProvider bridge={bridge}>
            <QueryProvider>
              <ApiProvider client={client} save={save}>
                <App />
              </ApiProvider>
            </QueryProvider>
          </BridgeProvider>
        </ErrorBoundary>
      </ThemeProvider>
    </StrictMode>,
  );
}
