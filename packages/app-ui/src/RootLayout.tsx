import { Outlet } from "@tanstack/react-router";

import { CommandNav } from "./components/CommandNav";
import { Header } from "./components/Header";
import { useElectronThemeSync } from "./hooks/useElectronThemeSync";
import { Tooltip } from "./ui";

export function RootLayout() {
  useElectronThemeSync();

  return (
    <Tooltip.Provider delay={300}>
      <div className="bg-surface-0 text-ink flex h-screen flex-col overflow-hidden">
        <Header />
        <div className="flex flex-1 overflow-hidden">
          <CommandNav />
          <Outlet />
        </div>
      </div>
    </Tooltip.Provider>
  );
}
