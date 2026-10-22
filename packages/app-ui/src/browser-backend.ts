import { createContext, useContext } from "react";

export interface BrowserBackend {
  readonly storage: "tab";
}

export const TAB_BACKEND: BrowserBackend = { storage: "tab" };

export const BrowserBackendContext = createContext<BrowserBackend | null>(null);

export function useBrowserBackend(): BrowserBackend | null {
  return useContext(BrowserBackendContext);
}
