import type { ApiClient } from "@neherlab/app-contracts/client";
import { useMemo } from "react";

import { ApiContext } from "./context";
import { useAppEvents } from "./events";

export function ApiProvider({ client, children }: { client: ApiClient; children: React.ReactNode }) {
  const value = useMemo(() => ({ client }), [client]);

  return (
    <ApiContext.Provider value={value}>
      <AppEvents />
      {children}
    </ApiContext.Provider>
  );
}

function AppEvents() {
  useAppEvents();

  return null;
}
