import { createContext, useContext } from "react";

import type { Host } from "./host";

export const HostContext = createContext<Host | null>(null);

export function useHost(): Host | null {
  return useContext(HostContext);
}
