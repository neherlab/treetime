import { useViewportSize } from "@mantine/hooks";
import { Outlet } from "@tanstack/react-router";

import { useNativeThemeSync } from "./hooks/useNativeThemeSync";
import { AppSidebar } from "./shell/AppSidebar";
import { CommandPalette } from "./shell/CommandPalette";
import { SiteHeader } from "./shell/SiteHeader";
import { useGlobalShortcuts } from "./shell/useGlobalShortcuts";
import { useYamlDrop } from "./shell/useYamlDrop";
import { WorkspaceDialog } from "./shell/WorkspaceDialog";
import { usePreferencesStore } from "./store/preferences";
import { SidebarInset, SidebarProvider } from "./ui/sidebar";
import { fitSidebarWidth, sidebarWidthOrDefault } from "./ui/sidebar-width";
import { TooltipProvider } from "./ui/tooltip";

export const MAIN_SCROLL_ID = "main-scroll";

export function RootLayout() {
  useNativeThemeSync();
  useGlobalShortcuts();
  const { getRootProps, getInputProps } = useYamlDrop();
  const { width: viewportWidth } = useViewportSize();
  const storedWidth = usePreferencesStore((state) => sidebarWidthOrDefault(state.preferences.sidebar_width));
  const sidebarWidth = fitSidebarWidth(storedWidth, viewportWidth);
  const setSidebarWidth = usePreferencesStore((state) => state.setSidebarWidth);

  return (
    <TooltipProvider delay={300}>
      <SidebarProvider
        width={sidebarWidth}
        onWidthChange={setSidebarWidth}
        {...getRootProps({ className: "h-svh flex-col overflow-hidden [--header-height:calc(--spacing(12))]" })}
      >
        <input {...getInputProps()} />
        <SiteHeader />
        <div className="flex min-h-0 flex-1">
          <AppSidebar />
          <SidebarInset
            id={MAIN_SCROLL_ID}
            data-scroll-restoration-id={MAIN_SCROLL_ID}
            className="@container min-h-0 overflow-y-auto overscroll-contain"
          >
            <Outlet />
          </SidebarInset>
        </div>
      </SidebarProvider>
      <CommandPalette />
      <WorkspaceDialog />
    </TooltipProvider>
  );
}
