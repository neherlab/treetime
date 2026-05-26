import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { App, BridgeProvider } from "@neherlab/treetime-ui";
import { createWebBridge } from "./bridge-web";

const bridge = createWebBridge();

createRoot(document.getElementById("root")!).render(
  <StrictMode>
    <BridgeProvider bridge={bridge}>
      <App />
    </BridgeProvider>
  </StrictMode>,
);
