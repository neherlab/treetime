import type { LocalFiles, TreeTimeBridge } from "@neherlab/app-contracts";
import { App, BridgeProvider, ErrorBoundary, QueryProvider, ThemeProvider } from "@neherlab/app-ui";
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";

import "./index.css";

declare global {
  interface Window {
    treetime: TreeTimeBridge;
    treetimeFiles: LocalFiles;
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
              <App localFiles={window.treetimeFiles} />
            </QueryProvider>
          </BridgeProvider>
        </ErrorBoundary>
      </ThemeProvider>
    </StrictMode>,
  );
}
