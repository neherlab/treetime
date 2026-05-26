import { createContext, useContext } from "react";
import type { TreeTimeBridge } from "@neherlab/app-contracts";

const BridgeContext = createContext<TreeTimeBridge | null>(null);

export function BridgeProvider({
  bridge,
  children,
}: {
  bridge: TreeTimeBridge;
  children: React.ReactNode;
}) {
  return <BridgeContext value={bridge}>{children}</BridgeContext>;
}

export function useBridge(): TreeTimeBridge {
  const bridge = useContext(BridgeContext);
  if (!bridge) {
    throw new Error("useBridge must be used within a BridgeProvider");
  }
  return bridge;
}
