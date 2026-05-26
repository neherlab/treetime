import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
<<<<<<< HEAD
import { App, BridgeProvider, QueryProvider } from "@neherlab/app-ui";
=======
import { App, BridgeProvider } from "@neherlab/app-ui";
>>>>>>> e3aa033b (feat(desktop): wire IPC handlers, add React renderer with Vite)
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
<<<<<<< HEAD
      <QueryProvider>
        <App />
      </QueryProvider>
=======
      <App />
>>>>>>> e3aa033b (feat(desktop): wire IPC handlers, add React renderer with Vite)
    </BridgeProvider>
  </StrictMode>,
);
