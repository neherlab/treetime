import { createContext, useContext } from "react";
import type { TreeTimeBridge } from "@neherlab/app-contracts";

export const BridgeContext = createContext<TreeTimeBridge | null>(null);

export function useBridge(): TreeTimeBridge {
  const bridge = useContext(BridgeContext);
  if (!bridge) {
    throw new Error("useBridge must be used within a BridgeProvider");
  }
  return bridge;
}
