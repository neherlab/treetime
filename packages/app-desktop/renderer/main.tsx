import { App, BridgeProvider, ErrorBoundary, QueryProvider, ThemeProvider } from "@neherlab/app-ui";
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";

import {
  createDesktopBridge,
  createLocalFiles,
  windowBackendConnection,
  type DesktopShell,
} from "../src/desktop-bridge";

import "./index.css";

declare global {
  interface Window {
    treetimeShell: DesktopShell;
  }
}

const bridge = createDesktopBridge(windowBackendConnection(window, window.treetimeShell), window.treetimeShell);

const localFiles = createLocalFiles(window.treetimeShell);

const root = document.getElementById("root");

if (root) {
  createRoot(root).render(
    <StrictMode>
      <ThemeProvider>
        <ErrorBoundary>
          <BridgeProvider bridge={bridge}>
            <QueryProvider>
              <App localFiles={localFiles} />
            </QueryProvider>
          </BridgeProvider>
        </ErrorBoundary>
      </ThemeProvider>
    </StrictMode>,
  );
}
