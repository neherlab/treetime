import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { App, BridgeProvider, QueryProvider } from "@neherlab/app-ui";
import type { TreeTimeBridge } from "@neherlab/app-contracts";

declare global {
  interface Window {
    treetime: TreeTimeBridge;
  }
}

const bridge = window.treetime;

createRoot(document.getElementById("root")!).render(
  <StrictMode>
    <BridgeProvider bridge={bridge}>
      <QueryProvider>
        <App />
      </QueryProvider>
    </BridgeProvider>
  </StrictMode>,
);
