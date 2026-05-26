import type { TreeTimeBridge } from "@neherlab/app-contracts";
import { BridgeContext } from "./BridgeContext";

<<<<<<< HEAD
export function BridgeProvider({ bridge, children }: { bridge: TreeTimeBridge; children: React.ReactNode }) {
=======
export function BridgeProvider({
  bridge,
  children,
}: {
  bridge: TreeTimeBridge;
  children: React.ReactNode;
}) {
>>>>>>> 8ba85e1a (fix(web): split BridgeContext for React Fast Refresh compatibility)
  return <BridgeContext value={bridge}>{children}</BridgeContext>;
}
