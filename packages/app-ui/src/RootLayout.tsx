import { Outlet } from "@tanstack/react-router";

import { useElectronThemeSync } from "./hooks/useElectronThemeSync";
import { CommandPalette } from "./shell/CommandPalette";
import { Sidebar } from "./shell/Sidebar";
import { TopBar } from "./shell/TopBar";
import { useGlobalShortcuts } from "./shell/useGlobalShortcuts";
import { useYamlDrop } from "./shell/useYamlDrop";
import { Tooltip } from "./ui";

export function RootLayout() {
  useElectronThemeSync();
  useGlobalShortcuts();
  useYamlDrop();

  return (
    <Tooltip.Provider delay={300}>
      <div className="bg-surface-0 text-ink isolate grid min-h-dvh grid-cols-[minmax(0,1fr)] grid-rows-[var(--spacing-bar)_1fr] text-sm">
        <TopBar />
        <div className="grid grid-cols-1 md:grid-cols-[18.75rem_minmax(0,1fr)]">
          <Sidebar />
          <main className="@container">
            <Outlet />
          </main>
        </div>
      </div>
      <CommandPalette />
    </Tooltip.Provider>
  );
}
