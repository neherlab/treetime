import type { ApiClient } from "@neherlab/app-contracts/client";
import { createContext, useContext } from "react";

export interface SaveActions {
  saveRunFile(id: string, path: string, name: string): Promise<boolean>;
  saveRunArchive(id: string, name: string): Promise<boolean>;
}

export interface ApiContextValue {
  client: ApiClient;
  save: SaveActions;
}

export const ApiContext = createContext<ApiContextValue | null>(null);

export function useApiContext(): ApiContextValue {
  const value = useContext(ApiContext);

  if (value === null) {
    throw new Error("useApiContext must be used within an ApiProvider");
  }

  return value;
}
