import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { App, BridgeProvider, QueryProvider } from "@neherlab/app-ui";
import { createMockBridge } from "./bridge-mock";
import "./index.css";

const bridge = createMockBridge();

const root = document.getElementById("root");
if (root) {
  createRoot(root).render(
    <StrictMode>
      <BridgeProvider bridge={bridge}>
        <QueryProvider>
          <App />
        </QueryProvider>
      </BridgeProvider>
    </StrictMode>,
  );
}
