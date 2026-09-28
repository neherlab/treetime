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
      <div className="bg-surface-0 text-ink relative grid h-screen grid-rows-[3rem_1fr] overflow-hidden text-sm">
        <TopBar />
        <div className="grid min-h-0 grid-cols-1 md:grid-cols-[18.75rem_1fr]">
          <Sidebar />
          <main className="relative min-w-0 overflow-auto">
            <Outlet />
          </main>
        </div>
      </div>
      <CommandPalette />
    </Tooltip.Provider>
  );
}
