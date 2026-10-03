import { mergeProps } from "@base-ui/react/merge-props";
import { useRender } from "@base-ui/react/use-render";
import { useMediaQuery } from "@mantine/hooks";
import { useHotkey } from "@tanstack/react-hotkeys";
import * as React from "react";
import PanelLeftIcon from "~icons/lucide/panel-left";

import { Button } from "./button";
import { cn } from "./cn";
import { Sheet, SheetContent, SheetDescription, SheetHeader, SheetTitle } from "./sheet";
import {
  clampSidebarWidth,
  SIDEBAR_WIDTH_DEFAULT,
  SIDEBAR_WIDTH_MAX,
  SIDEBAR_WIDTH_MIN,
  sidebarWidthForKey,
} from "./sidebar-width";

const MOBILE_QUERY = "(max-width: 767px)";

const SIDEBAR_HOTKEY = "Mod+B";

declare module "react" {
  interface CSSProperties {
    "--sidebar-width"?: string;
  }
}

const SidebarContext = React.createContext<SidebarContextProps | undefined>(undefined);

function SidebarProvider({
  width,
  onWidthChange: setWidth,
  className,
  style,
  children,
  ...props
}: React.ComponentProps<"div"> & { width: number; onWidthChange: (width: number) => void }) {
  const isMobile = useMediaQuery(MOBILE_QUERY);
  const [open, setOpen] = React.useState(true);
  const [openMobile, setOpenMobile] = React.useState(false);
  const [resizing, setResizing] = React.useState(false);

  const toggleSidebar = React.useCallback(() => {
    if (isMobile) {
      setOpenMobile((value) => !value);
    } else {
      setOpen((value) => !value);
    }
  }, [isMobile]);

  useHotkey(SIDEBAR_HOTKEY, toggleSidebar);

  const contextValue = React.useMemo<SidebarContextProps>(
    () => ({ open, isMobile, openMobile, setOpenMobile, toggleSidebar, width, setWidth, resizing, setResizing }),
    [open, isMobile, openMobile, toggleSidebar, width, setWidth, resizing],
  );

  return (
    <SidebarContext value={contextValue}>
      <div
        data-slot="sidebar-wrapper"
        className={cn("group/sidebar-wrapper flex min-h-svh w-full", className)}
        // oxlint-disable-next-line react/forbid-dom-props -- the dragged sidebar width is a runtime value that no Tailwind class can hold
        style={{ ...style, "--sidebar-width": `${width}px` }}
        {...props}
      >
        {children}
      </div>
    </SidebarContext>
  );
}

function Sidebar({ className, children }: { className?: string; children: React.ReactNode }) {
  const { isMobile, open, openMobile, setOpenMobile, resizing } = useSidebar();

  if (isMobile) {
    return (
      <Sheet open={openMobile} onOpenChange={setOpenMobile}>
        <SheetContent
          data-slot="sidebar"
          data-mobile="true"
          className="bg-sidebar text-sidebar-foreground w-(--sidebar-width) p-0 [--sidebar-width:18rem] [&>button]:hidden"
          side="left"
        >
          <SheetHeader className="sr-only">
            <SheetTitle>Sidebar</SheetTitle>
            <SheetDescription>Displays the mobile sidebar.</SheetDescription>
          </SheetHeader>
          <div className="flex h-full w-full flex-col">{children}</div>
        </SheetContent>
      </Sheet>
    );
  }

  return (
    <div
      className="group peer text-sidebar-foreground hidden md:block"
      data-state={open ? "expanded" : "collapsed"}
      data-resizing={resizing || undefined}
      data-slot="sidebar"
    >
      <div
        data-slot="sidebar-gap"
        className="relative w-(--sidebar-width) bg-transparent transition-[width] duration-200 ease-linear group-data-resizing:transition-none group-data-[state=collapsed]:w-0"
      />
      <div
        data-slot="sidebar-container"
        className={cn(
          "fixed inset-y-0 left-0 z-10 hidden h-svh w-(--sidebar-width) border-r transition-[left] duration-200 ease-linear group-data-resizing:transition-none group-data-[state=collapsed]:left-[calc(var(--sidebar-width)*-1)] md:flex",
          className,
        )}
      >
        <div data-slot="sidebar-inner" className="bg-sidebar flex size-full flex-col">
          {children}
        </div>
        <SidebarResizeHandle />
      </div>
    </div>
  );
}

function SidebarResizeHandle() {
  const { width, setWidth, resizing, setResizing } = useSidebar();

  const onPointerDown = React.useCallback(
    (event: React.PointerEvent<HTMLDivElement>) => {
      if (event.button !== 0) {
        return;
      }

      event.preventDefault();
      event.currentTarget.setPointerCapture(event.pointerId);
      setResizing(true);
    },
    [setResizing],
  );

  const onPointerMove = React.useCallback(
    (event: React.PointerEvent<HTMLDivElement>) => {
      const container = event.currentTarget.parentElement;

      if (container !== null && event.currentTarget.hasPointerCapture(event.pointerId)) {
        setWidth(clampSidebarWidth(event.clientX - container.getBoundingClientRect().left));
      }
    },
    [setWidth],
  );

  const onLostPointerCapture = React.useCallback(() => setResizing(false), [setResizing]);

  const onKeyDown = React.useCallback(
    (event: React.KeyboardEvent<HTMLDivElement>) => {
      const next = sidebarWidthForKey(width, event.key);

      if (next !== undefined) {
        event.preventDefault();
        setWidth(next);
      }
    },
    [setWidth, width],
  );

  const onDoubleClick = React.useCallback(() => setWidth(SIDEBAR_WIDTH_DEFAULT), [setWidth]);

  return (
    <div
      data-slot="sidebar-resize-handle"
      data-resizing={resizing || undefined}
      // oxlint-disable-next-line jsx-a11y/prefer-tag-over-role -- a focusable separator is the WAI-ARIA window splitter widget, and the linter treats hr as non-interactive
      role="separator"
      aria-orientation="vertical"
      aria-label="Resize sidebar"
      aria-valuenow={width}
      aria-valuemin={SIDEBAR_WIDTH_MIN}
      aria-valuemax={SIDEBAR_WIDTH_MAX}
      tabIndex={0}
      title="Drag to resize, double-click to reset"
      onPointerDown={onPointerDown}
      onPointerMove={onPointerMove}
      onLostPointerCapture={onLostPointerCapture}
      onKeyDown={onKeyDown}
      onDoubleClick={onDoubleClick}
      className="after:bg-sidebar-ring absolute inset-y-0 -right-1 z-20 m-0 h-auto w-2 cursor-col-resize touch-none border-0 outline-hidden group-data-[state=collapsed]:hidden after:absolute after:inset-y-0 after:left-1/2 after:w-0.5 after:-translate-x-1/2 after:opacity-0 after:transition-opacity hover:after:opacity-100 focus-visible:after:opacity-100 data-resizing:after:opacity-100"
    />
  );
}

function SidebarTrigger({ className, ...props }: Omit<React.ComponentProps<typeof Button>, "onClick">) {
  const { toggleSidebar } = useSidebar();

  return (
    <Button
      data-slot="sidebar-trigger"
      variant="ghost"
      size="icon-sm"
      className={className}
      onClick={toggleSidebar}
      {...props}
    >
      <PanelLeftIcon aria-hidden />
      <span className="sr-only">Toggle Sidebar</span>
    </Button>
  );
}

function SidebarInset({ className, ...props }: React.ComponentProps<"main">) {
  return (
    <main
      data-slot="sidebar-inset"
      className={cn("bg-background relative flex w-full flex-1 flex-col", className)}
      {...props}
    />
  );
}

function SidebarHeader({ className, ...props }: React.ComponentProps<"div">) {
  return <div data-slot="sidebar-header" className={cn("flex flex-col gap-2 p-2", className)} {...props} />;
}

function SidebarFooter({ className, ...props }: React.ComponentProps<"div">) {
  return <div data-slot="sidebar-footer" className={cn("flex flex-col gap-2 p-2", className)} {...props} />;
}

function SidebarContent({ className, ...props }: React.ComponentProps<"div">) {
  return (
    <div
      data-slot="sidebar-content"
      className={cn("no-scrollbar flex min-h-0 flex-1 flex-col gap-2 overflow-auto", className)}
      {...props}
    />
  );
}

function SidebarGroup({ className, ...props }: React.ComponentProps<"div">) {
  return (
    <div data-slot="sidebar-group" className={cn("relative flex w-full min-w-0 flex-col p-2", className)} {...props} />
  );
}

function SidebarGroupLabel({
  className,
  render,
  ...props
}: useRender.ComponentProps<"div"> & React.ComponentProps<"div">) {
  return useRender({
    defaultTagName: "div",
    props: mergeProps<"div">(
      {
        className: cn(
          "text-sidebar-foreground/70 ring-sidebar-ring flex h-8 shrink-0 items-center rounded-md px-2 text-xs font-bold outline-hidden focus-visible:ring-2 [&>svg]:size-4 [&>svg]:shrink-0",
          className,
        ),
      },
      props,
    ),
    render,
    state: { slot: "sidebar-group-label" },
  });
}

function SidebarMenu({ className, ...props }: React.ComponentProps<"ul">) {
  return <ul data-slot="sidebar-menu" className={cn("flex w-full min-w-0 flex-col gap-1", className)} {...props} />;
}

function SidebarMenuItem({ className, ...props }: React.ComponentProps<"li">) {
  return <li data-slot="sidebar-menu-item" className={cn("group/menu-item relative", className)} {...props} />;
}

function SidebarMenuButton({
  render,
  isActive = false,
  className,
  ...props
}: useRender.ComponentProps<"button"> & React.ComponentProps<"button"> & { isActive?: boolean }) {
  return useRender({
    defaultTagName: "button",
    props: mergeProps<"button">(
      {
        className: cn(
          "peer/menu-button group/menu-button ring-sidebar-ring hover:bg-sidebar-accent hover:text-sidebar-accent-foreground active:bg-sidebar-accent active:text-sidebar-accent-foreground data-active:bg-sidebar-accent data-active:text-sidebar-accent-foreground flex h-8 w-full items-center gap-2 overflow-hidden rounded-md p-2 text-left text-sm outline-hidden transition-[width,height,padding] focus-visible:ring-2 disabled:pointer-events-none disabled:opacity-50 aria-disabled:pointer-events-none aria-disabled:opacity-50 [&_svg]:size-4 [&_svg]:shrink-0 [&>span:last-child]:truncate",
          className,
        ),
      },
      props,
    ),
    render,
    state: { slot: "sidebar-menu-button", active: isActive },
  });
}

interface SidebarContextProps {
  open: boolean;
  isMobile: boolean;
  openMobile: boolean;
  setOpenMobile: (open: boolean) => void;
  toggleSidebar: () => void;
  width: number;
  setWidth: (width: number) => void;
  resizing: boolean;
  setResizing: (resizing: boolean) => void;
}

function useSidebar(): SidebarContextProps {
  const context = React.use(SidebarContext);

  if (context === undefined) {
    throw new Error("useSidebar must be used within a SidebarProvider.");
  }

  return context;
}

export {
  Sidebar,
  SidebarContent,
  SidebarFooter,
  SidebarGroup,
  SidebarGroupLabel,
  SidebarHeader,
  SidebarInset,
  SidebarMenu,
  SidebarMenuButton,
  SidebarMenuItem,
  SidebarProvider,
  SidebarTrigger,
};
