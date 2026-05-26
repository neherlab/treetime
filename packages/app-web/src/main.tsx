import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { App, BridgeProvider, QueryProvider, ThemeProvider } from "@neherlab/app-ui";
import { createWebBridge } from "./bridge-web";
import "./index.css";

const bridge = createWebBridge();

const root = document.getElementById("root");
if (root) {
  createRoot(root).render(
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
}
