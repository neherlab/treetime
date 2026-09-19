import { App, BridgeProvider, ErrorBoundary, QueryProvider, ThemeProvider } from "@neherlab/app-ui";
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";

import { createWebBridge } from "./bridge-web";

import "./index.css";

const bridge = createWebBridge();

const root = document.getElementById("root");
if (root) {
  createRoot(root).render(
    <StrictMode>
      <ThemeProvider>
        <ErrorBoundary>
          <BridgeProvider bridge={bridge}>
            <QueryProvider>
              <App />
            </QueryProvider>
          </BridgeProvider>
        </ErrorBoundary>
      </ThemeProvider>
    </StrictMode>,
  );
}
