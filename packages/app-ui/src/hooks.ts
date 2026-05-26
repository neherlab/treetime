import type { VersionInfo } from "@neherlab/app-contracts";
import { useQuery } from "@tanstack/react-query";
import { useBridge } from "./BridgeContext";

export function useVersion() {
  const bridge = useBridge();
  return useQuery<VersionInfo>({
    queryKey: ["version"],
    queryFn: () => bridge.version(),
    staleTime: Infinity,
  });
}
