import { Dialog } from "@base-ui/react/dialog";
import { useReducedMotion } from "@mantine/hooks";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "@neherlab/design/ui/collapsible";
import { CHANGE_PANEL_LAYOUT, TOGGLE_SIDEBAR } from "auspice/src/actions/types";
import {
  useCallback,
  useEffect,
  useId,
  useRef,
  useState,
  type CSSProperties,
  type ReactNode,
  type RefObject,
} from "react";

import { AuspicePanels } from "./AuspicePanels";
import { AuspiceProviders } from "./AuspiceProviders";
import { HideControlsButton, SidebarControls } from "./AuspiceSidebar";
import type { FigureRuns } from "./downloads";
import { isSidebarDocked, SIDEBAR_WIDTH, startPanelLayout } from "./layout";
import { useNotificationToasts } from "./notifications";
import { observeHeight, scrollAncestor } from "./scrollPane";
import type { AuspiceState } from "./state";
import type { AuspiceStore } from "./store";
import { useAuspiceSelector } from "./store-hooks";
import { ViewerToolbar } from "./ViewerToolbar";

const SLIDE_MS = 300;

const STARTED_STORES = new WeakSet<AuspiceStore>();

const SIDEBAR_CLASS =
  "bg-muted absolute inset-y-0 left-0 z-20 w-(--auspice-sidebar-width) overflow-x-hidden overflow-y-auto border-r transition-[left] duration-300 ease-out motion-reduce:transition-none data-ending-style:-left-(--auspice-sidebar-width) data-starting-style:-left-(--auspice-sidebar-width)";

const VIEWER_STYLE: CSSProperties = { "--auspice-sidebar-width": `${SIDEBAR_WIDTH}px` };

declare module "react" {
  interface CSSProperties {
    "--auspice-sidebar-width"?: string;
  }
}

export function AuspiceTree({
  store,
  runs,
  version,
  entropy = true,
}: {
  store: AuspiceStore;
  runs: FigureRuns;
  version: string | undefined;
  entropy?: boolean | undefined;
}) {
  useNotificationToasts(store);

  return (
    <AuspiceProviders store={store}>
      <Viewer store={store} runs={runs} version={version} entropy={entropy} />
    </AuspiceProviders>
  );
}

function Viewer({
  store,
  runs,
  version,
  entropy,
}: {
  store: AuspiceStore;
  runs: FigureRuns;
  version: string | undefined;
  entropy: boolean;
}) {
  const viewerRef = useRef<HTMLDivElement | null>(null);
  const [measureViewer, viewerWidth] = useViewerWidth(store, viewerRef);
  const visibleHeight = useVisibleHeight(viewerRef);
  const open = useAuspiceSelector(store, selectSidebarOpen);
  const docked = isSidebarDocked(viewerWidth);
  const { slideFrom, setOpen } = useSidebarToggle(store, docked);
  const sidebarId = useId();
  const sizedOpen = docked && (slideFrom ?? open);

  const panels = (
    <AuspicePanels
      store={store}
      panelsWidth={viewerWidth - (sizedOpen ? SIDEBAR_WIDTH : 0)}
      visibleHeight={visibleHeight}
      entropy={entropy}
      runs={runs}
      version={version}
    />
  );

  const optionsRef = useRef<HTMLButtonElement>(null);

  return (
    <div className="grid min-w-0 gap-2">
      <ViewerToolbar
        store={store}
        sidebarId={sidebarId}
        sidebarOpen={open}
        onSidebarOpen={setOpen}
        optionsRef={optionsRef}
      />
      <div
        ref={measureViewer}
        // oxlint-disable-next-line react/forbid-dom-props -- the sidebar width comes from Auspice's layout constants, which no Tailwind class can hold
        style={VIEWER_STYLE}
        className="bg-card text-foreground relative min-w-0 overflow-hidden rounded-lg border scheme-light"
      >
        {docked && (
          <DockedSidebar
            store={store}
            open={open}
            onOpenChange={setOpen}
            sidebarId={sidebarId}
            entropy={entropy}
            optionsRef={optionsRef}
          >
            {panels}
          </DockedSidebar>
        )}
        {viewerWidth > 0 && !docked && (
          <OverlaySidebar
            store={store}
            open={open}
            onOpenChange={setOpen}
            sidebarId={sidebarId}
            entropy={entropy}
            optionsRef={optionsRef}
            viewerRef={viewerRef}
          >
            {panels}
          </OverlaySidebar>
        )}
      </div>
    </div>
  );
}

function DockedSidebar({
  store,
  open,
  onOpenChange,
  sidebarId,
  entropy,
  optionsRef,
  children,
}: SidebarProps & { children: ReactNode }) {
  const sidebarRef = useRef<HTMLDivElement>(null);

  const changeOpen = useCallback(
    (value: boolean) => {
      const hadFocus = sidebarRef.current?.contains(document.activeElement) === true;

      onOpenChange(value);

      if (!value && hadFocus) {
        optionsRef.current?.focus();
      }
    },
    [onOpenChange, optionsRef],
  );

  return (
    <Collapsible open={open} onOpenChange={changeOpen}>
      <CollapsibleContent keepMounted id={sidebarId} ref={sidebarRef} className={SIDEBAR_CLASS}>
        <CollapsibleTrigger render={<HideControlsButton />} />
        <SidebarControls store={store} entropy={entropy} />
      </CollapsibleContent>
      <div
        data-open={open ? "" : undefined}
        className="relative isolate transition-[margin-left] duration-300 ease-out data-open:ml-(--auspice-sidebar-width) motion-reduce:transition-none"
      >
        {children}
      </div>
    </Collapsible>
  );
}

function OverlaySidebar({
  store,
  open,
  onOpenChange,
  sidebarId,
  entropy,
  optionsRef,
  viewerRef,
  children,
}: SidebarProps & { viewerRef: RefObject<HTMLDivElement | null>; children: ReactNode }) {
  return (
    <Dialog.Root open={open} onOpenChange={onOpenChange} modal="trap-focus">
      <div className="relative isolate">{children}</div>
      <Dialog.Portal container={viewerRef}>
        <Dialog.Backdrop className="absolute inset-y-0 right-0 left-(--auspice-sidebar-width) z-10 bg-black/50 transition-opacity duration-300 ease-out data-ending-style:opacity-0 data-starting-style:opacity-0 motion-reduce:transition-none" />
        <Dialog.Popup id={sidebarId} aria-label="Auspice controls" className={SIDEBAR_CLASS} finalFocus={optionsRef}>
          <Dialog.Close render={<HideControlsButton />} />
          <SidebarControls store={store} entropy={entropy} />
        </Dialog.Popup>
      </Dialog.Portal>
    </Dialog.Root>
  );
}

interface SidebarProps {
  store: AuspiceStore;
  open: boolean;
  onOpenChange: (open: boolean) => void;
  sidebarId: string;
  entropy: boolean;
  optionsRef: RefObject<HTMLButtonElement | null>;
}

function useViewerWidth(
  store: AuspiceStore,
  viewerRef: RefObject<HTMLDivElement | null>,
): [(viewer: HTMLDivElement | null) => () => void, number] {
  const [width, setWidth] = useState(0);

  const measure = useCallback(
    (viewer: HTMLDivElement | null) => {
      viewerRef.current = viewer;

      let started = false;

      const observer = new ResizeObserver(([entry]) => {
        const measured = Math.floor(entry?.contentRect.width ?? 0);

        if (!started && measured > 0) {
          started = true;
          applyStartLayout(store, measured);
        }

        setWidth(measured);
      });

      if (viewer !== null) {
        observer.observe(viewer);
      }

      return () => {
        observer.disconnect();
        viewerRef.current = null;
      };
    },
    [store, viewerRef],
  );

  return [measure, width];
}

function applyStartLayout(store: AuspiceStore, viewerWidth: number): void {
  if (STARTED_STORES.has(store)) {
    return;
  }

  STARTED_STORES.add(store);

  store.dispatch({ type: TOGGLE_SIDEBAR, value: false });

  if (store.getState().controls.canTogglePanelLayout) {
    store.dispatch({ type: CHANGE_PANEL_LAYOUT, data: startPanelLayout(viewerWidth) });
  }
}

function useVisibleHeight(viewerRef: RefObject<HTMLDivElement | null>): number {
  const [height, setHeight] = useState(() => window.innerHeight);

  useEffect(() => observeHeight(scrollAncestor(viewerRef.current), setHeight), [viewerRef]);

  return height;
}

function useSidebarToggle(store: AuspiceStore, docked: boolean): SidebarToggle {
  const reducedMotion = useReducedMotion();
  const [slideFrom, setSlideFrom] = useState<boolean | undefined>(undefined);
  const [timer] = useState(() => new SlideTimer());

  useEffect(() => () => timer.stop(), [timer]);

  const setOpen = useCallback(
    (value: boolean) => {
      const current = store.getState().controls.sidebarOpen;

      if (value === current) {
        return;
      }

      if (docked && !reducedMotion) {
        setSlideFrom((from) => from ?? current);
        timer.start(() => setSlideFrom(undefined));
      }

      store.dispatch({ type: TOGGLE_SIDEBAR, value });
    },
    [docked, reducedMotion, store, timer],
  );

  return { slideFrom: docked ? slideFrom : undefined, setOpen };
}

interface SidebarToggle {
  slideFrom: boolean | undefined;
  setOpen: (open: boolean) => void;
}

class SlideTimer {
  #timeout: number | undefined;

  start(onEnd: () => void): void {
    this.stop();
    this.#timeout = window.setTimeout(onEnd, SLIDE_MS);
  }

  stop(): void {
    window.clearTimeout(this.#timeout);
    this.#timeout = undefined;
  }
}

function selectSidebarOpen(state: AuspiceState): boolean {
  return state.controls.sidebarOpen;
}
