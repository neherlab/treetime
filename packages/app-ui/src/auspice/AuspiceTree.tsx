import { Dialog } from "@base-ui/react/dialog";
import { useReducedMotion } from "@mantine/hooks";
import { CHANGE_PANEL_LAYOUT, TOGGLE_SIDEBAR } from "auspice/src/actions/types";
import {
  useCallback,
  useEffect,
  useId,
  useRef,
  useState,
  useSyncExternalStore,
  type CSSProperties,
  type ReactNode,
  type RefObject,
} from "react";
import { I18nextProvider } from "react-i18next";
import { Provider } from "react-redux";
import { ThemeProvider } from "styled-components";

import { MAIN_SCROLL_ID } from "../RootLayout";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "../ui/collapsible";
import { AuspicePanels } from "./AuspicePanels";
import { HideControlsButton, ShowControlsTab, SidebarControls } from "./AuspiceSidebar";
import type { FigureRuns } from "./downloads";
import { AUSPICE_I18N } from "./i18n";
import { isSidebarDocked, SIDEBAR_WIDTH, startLayout } from "./layout";
import { useNotificationToasts } from "./notifications";
import type { AuspiceState } from "./state";
import type { AuspiceStore } from "./store";
import { useAuspiceSelector } from "./store-hooks";

const SIDEBAR_THEME = {
  background: "var(--muted)",
  color: "var(--foreground)",
  "font-family": "Lato, Helvetica Neue, Helvetica, sans-serif",
  sidebarBoxShadow: "rgba(0, 0, 0, 0.15)",
  selectedColor: "var(--primary)",
  unselectedColor: "var(--muted-foreground)",
  alternateBackground: "var(--input)",
};

const SLIDE_MS = 300;

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
    <I18nextProvider i18n={AUSPICE_I18N}>
      <ThemeProvider theme={SIDEBAR_THEME}>
        <Provider store={store}>
          <Viewer store={store} runs={runs} version={version} entropy={entropy} />
        </Provider>
      </ThemeProvider>
    </I18nextProvider>
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
  const visibleHeight = useVisibleHeight();
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

  return (
    <div
      ref={measureViewer}
      // oxlint-disable-next-line react/forbid-dom-props -- the sidebar width comes from Auspice's layout constants, which no Tailwind class can hold
      style={VIEWER_STYLE}
      className="light-scope relative min-w-0 overflow-hidden rounded-lg border bg-white"
    >
      {docked && (
        <DockedSidebar store={store} open={open} onOpenChange={setOpen} sidebarId={sidebarId} entropy={entropy}>
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
          viewerRef={viewerRef}
        >
          {panels}
        </OverlaySidebar>
      )}
    </div>
  );
}

function DockedSidebar({
  store,
  open,
  onOpenChange,
  sidebarId,
  entropy,
  children,
}: SidebarProps & { children: ReactNode }) {
  const sidebarRef = useRef<HTMLDivElement>(null);
  const focusTab = useRef(false);

  const changeOpen = useCallback(
    (value: boolean) => {
      focusTab.current = !value && sidebarRef.current?.contains(document.activeElement) === true;
      onOpenChange(value);
    },
    [onOpenChange],
  );

  const tabRef = useCallback((tab: HTMLButtonElement | null) => {
    if (tab !== null && focusTab.current) {
      focusTab.current = false;
      tab.focus();
    }
  }, []);

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
        {!open && <CollapsibleTrigger render={<ShowControlsTab ref={tabRef} />} />}
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
  viewerRef,
  children,
}: SidebarProps & { viewerRef: RefObject<HTMLDivElement | null>; children: ReactNode }) {
  return (
    <Dialog.Root open={open} onOpenChange={onOpenChange} modal="trap-focus">
      <div className="relative isolate">
        <Dialog.Trigger render={<ShowControlsTab />} />
        {children}
      </div>
      <Dialog.Portal container={viewerRef}>
        <Dialog.Backdrop className="absolute inset-y-0 right-0 left-(--auspice-sidebar-width) z-10 bg-black/50 transition-opacity duration-300 ease-out data-ending-style:opacity-0 data-starting-style:opacity-0 motion-reduce:transition-none" />
        <Dialog.Popup id={sidebarId} aria-label="Auspice controls" className={SIDEBAR_CLASS}>
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
  const { sidebarOpen, panelLayout } = startLayout(viewerWidth);

  store.dispatch({ type: TOGGLE_SIDEBAR, value: sidebarOpen });

  if (store.getState().controls.canTogglePanelLayout) {
    store.dispatch({ type: CHANGE_PANEL_LAYOUT, data: panelLayout });
  }
}

function useVisibleHeight(): number {
  return useSyncExternalStore(subscribeVisibleHeight, visibleHeight, visibleHeight);
}

function subscribeVisibleHeight(onChange: () => void): () => void {
  const scroll = document.getElementById(MAIN_SCROLL_ID);

  if (scroll === null) {
    window.addEventListener("resize", onChange);

    return () => window.removeEventListener("resize", onChange);
  }

  const observer = new ResizeObserver(onChange);

  observer.observe(scroll);

  return () => observer.disconnect();
}

function visibleHeight(): number {
  return document.getElementById(MAIN_SCROLL_ID)?.clientHeight ?? window.innerHeight;
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
