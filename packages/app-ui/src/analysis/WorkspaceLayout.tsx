import { useMediaQuery } from "@mantine/hooks";
import { Button } from "@neherlab/design/ui/button";
import { type PanelSize, ResizableHandle, ResizablePanel, ResizablePanelGroup } from "@neherlab/design/ui/resizable";
import { Sheet, SheetContent, SheetDescription, SheetHeader, SheetTitle } from "@neherlab/design/ui/sheet";
import { useCallback, useRef, useState, type ReactNode } from "react";
import SlidersHorizontal from "~icons/lucide/sliders-horizontal";

import { usePreferencesStore } from "../store/preferences";

export const WORKSPACE_MAIN_ID = "workspace-main";

const PHONE_QUERY = "(max-width: 1023px)";

const DEFAULT_OPTIONS_WIDTH = "26rem";

export function WorkspaceLayout({
  title,
  description,
  actions,
  panel,
  footer,
  aboveMain,
  main,
}: {
  title: ReactNode;
  description: ReactNode;
  actions: ReactNode;
  panel: ReactNode;
  footer: ReactNode;
  aboveMain?: ReactNode;
  main: ReactNode;
}) {
  const phone = useMediaQuery(PHONE_QUERY, false, { getInitialValueInEffect: false });
  const [sheetOpen, setSheetOpen] = useState(false);
  const openSheet = useCallback(() => setSheetOpen(true), []);

  const mainArea = (
    <div
      id={WORKSPACE_MAIN_ID}
      data-scroll-restoration-id={WORKSPACE_MAIN_ID}
      className="@container h-full min-h-0 overflow-y-auto overscroll-contain"
    >
      {aboveMain}
      {main}
    </div>
  );

  const header = (
    <header className="flex flex-wrap items-center gap-x-4 gap-y-1 border-b px-4 py-2.5">
      <div className="grid min-w-0 gap-0.5">
        <h1 className="font-heading truncate text-lg leading-tight font-bold">{title}</h1>
        {description}
      </div>
      <div className="ml-auto flex flex-wrap items-center gap-1.5">
        {phone && (
          <Button type="button" variant="outline" size="sm" onClick={openSheet}>
            <SlidersHorizontal aria-hidden />
            Options
          </Button>
        )}
        {actions}
      </div>
    </header>
  );

  if (phone) {
    return (
      <div className="flex h-full min-h-0 flex-col">
        {header}
        <div className="min-h-0 flex-1">{mainArea}</div>
        <div className="bg-background border-t px-4 py-2.5">{footer}</div>
        <Sheet open={sheetOpen} onOpenChange={setSheetOpen}>
          <SheetContent side="left" className="w-[min(100%,28rem)] gap-0 p-0 data-[side=left]:sm:max-w-md">
            <SheetHeader className="border-b">
              <SheetTitle>Options</SheetTitle>
              <SheetDescription>{title}</SheetDescription>
            </SheetHeader>
            <div className="@container min-h-0 flex-1 overflow-y-auto overscroll-contain">{panel}</div>
          </SheetContent>
        </Sheet>
      </div>
    );
  }

  return (
    <div className="flex h-full min-h-0 flex-col">
      {header}
      <OptionsSplit panel={panel} footer={footer} main={mainArea} />
    </div>
  );
}

function OptionsSplit({ panel, footer, main }: { panel: ReactNode; footer: ReactNode; main: ReactNode }) {
  const [savedWidth] = useState(() => usePreferencesStore.getState().preferences.options_width);
  const width = useRef<number | undefined>(undefined);

  const onResize = useCallback((size: PanelSize) => {
    width.current = Math.round(size.inPixels);
  }, []);

  const onLayoutChanged = useCallback(() => {
    if (width.current !== undefined && width.current !== usePreferencesStore.getState().preferences.options_width) {
      usePreferencesStore.getState().setOptionsWidth(width.current);
    }
  }, []);

  return (
    <ResizablePanelGroup orientation="horizontal" className="min-h-0 flex-1" onLayoutChanged={onLayoutChanged}>
      <ResizablePanel
        id="options"
        defaultSize={savedWidth ?? DEFAULT_OPTIONS_WIDTH}
        minSize="20rem"
        maxSize="50%"
        groupResizeBehavior="preserve-pixel-size"
        onResize={onResize}
      >
        <aside aria-label="Options" className="flex h-full min-h-0 flex-col">
          <div className="@container min-h-0 flex-1 overflow-y-auto overscroll-contain">{panel}</div>
          <div className="border-t px-4 py-3">{footer}</div>
        </aside>
      </ResizablePanel>
      <ResizableHandle aria-label="Width of the options" />
      <ResizablePanel id="main" minSize="20rem">
        {main}
      </ResizablePanel>
    </ResizablePanelGroup>
  );
}
