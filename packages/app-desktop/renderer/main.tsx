import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { App, BridgeProvider, QueryProvider, ThemeProvider } from "@neherlab/app-ui";
import type { TreeTimeBridge } from "@neherlab/app-contracts";
import "./index.css";

declare global {
  interface Window {
    treetime: TreeTimeBridge;
  }
}

const bridge = window.treetime;

createRoot(document.getElementById("root")!).render(
  <StrictMode>
    <ThemeProvider>
      <BridgeProvider bridge={bridge}>
        <QueryProvider>
          <App />
        </QueryProvider>
      </BridgeProvider>
    </ThemeProvider>
  </StrictMode>,
);
