import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
<<<<<<< HEAD
<<<<<<< HEAD
import { App, BridgeProvider, QueryProvider } from "@neherlab/app-ui";
=======
import { App, BridgeProvider } from "@neherlab/app-ui";
>>>>>>> e3aa033b (feat(desktop): wire IPC handlers, add React renderer with Vite)
=======
import { App, BridgeProvider, QueryProvider, ThemeProvider } from "@neherlab/app-ui";
>>>>>>> 0571c0c9 (feat(ui): add automatic OS theme detection with next-themes)
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
<<<<<<< HEAD
    <BridgeProvider bridge={bridge}>
<<<<<<< HEAD
      <QueryProvider>
        <App />
      </QueryProvider>
=======
      <App />
>>>>>>> e3aa033b (feat(desktop): wire IPC handlers, add React renderer with Vite)
    </BridgeProvider>
=======
    <ThemeProvider>
      <BridgeProvider bridge={bridge}>
        <QueryProvider>
          <App />
        </QueryProvider>
      </BridgeProvider>
    </ThemeProvider>
>>>>>>> 0571c0c9 (feat(ui): add automatic OS theme detection with next-themes)
  </StrictMode>,
);
