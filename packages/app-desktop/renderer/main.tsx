import type { TreeTimeBridge } from "@neherlab/app-contracts";
import { App, BridgeProvider, ErrorBoundary, QueryProvider, ThemeProvider } from "@neherlab/app-ui";
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";

import "./index.css";

declare global {
  interface Window {
    treetime: TreeTimeBridge;
  }
}

const bridge = window.treetime;

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
