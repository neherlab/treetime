import type { ApiClient } from "@neherlab/app-contracts/client";
import { createContext, useContext } from "react";

export interface ApiContextValue {
  client: ApiClient;
}

export const ApiContext = createContext<ApiContextValue | null>(null);

export function useApiContext(): ApiContextValue {
  const value = useContext(ApiContext);

  if (value === null) {
    throw new Error("useApiContext must be used within an ApiProvider");
  }

  return value;
}
