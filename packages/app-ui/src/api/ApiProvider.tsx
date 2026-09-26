import type { ApiClient } from "@neherlab/app-contracts/client";
import { useMemo } from "react";

import { ApiContext, type SaveActions } from "./context";
import { useAppEvents } from "./events";

export function ApiProvider({
  client,
  save,
  children,
}: {
  client: ApiClient;
  save: SaveActions;
  children: React.ReactNode;
}) {
  const value = useMemo(() => ({ client, save }), [client, save]);

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
