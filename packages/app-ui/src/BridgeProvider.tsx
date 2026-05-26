import type { TreeTimeBridge } from "@neherlab/app-contracts";
import { BridgeContext } from "./BridgeContext";

export function BridgeProvider({ bridge, children }: { bridge: TreeTimeBridge; children: React.ReactNode }) {
  return <BridgeContext value={bridge}>{children}</BridgeContext>;
}
