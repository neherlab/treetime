import { useElementSize } from "@mantine/hooks";
import { errorMessage } from "@neherlab/app-contracts";
import Entropy from "auspice/src/components/entropy";
import FiltersSummary from "auspice/src/components/info/filtersSummary";
import Tree from "auspice/src/components/tree";
import {
  lazy,
  Suspense,
  useEffect,
  useLayoutEffect,
  useMemo,
  useRef,
  useState,
  type CSSProperties,
  type RefCallback,
  type RefObject,
} from "react";
import { ErrorBoundary, type FallbackProps } from "react-error-boundary";

import { Alert, AlertDescription, AlertTitle } from "../ui/alert";
import type { FigureRuns } from "./downloads";
import {
  CARD_TITLE_HEIGHT,
  CARD_TITLE_HEIGHT_TWO_TREES,
  panelSizes,
  type DrawingSize,
  type PanelSizes,
} from "./layout";
import type { AuspicePanelLayout, AuspiceState } from "./state";
import type { AuspiceStore } from "./store";
import { useAuspiceSelector } from "./store-hooks";
import { TreeDownloads } from "./TreeDownloads";
import { canDrawMaps } from "./webgl";

const AuspiceMap = lazy(async () => import("auspice/src/components/map/map"));

const ZOOM_BUTTONS = ":scope > div:last-child > div + div > div:last-child";

const ZOOM_BUTTONS_GAP = 4;

const HEADER_ITEM_GAP = 16;

const HEADER_MIN_SUMMARY_WIDTH = 120;

declare module "react" {
  interface CSSProperties {
    "--auspice-card-title-height"?: string;
    "--auspice-header-right"?: string;
  }
}

export function AuspicePanels({
  store,
  panelsWidth,
  visibleHeight,
  entropy,
  runs,
  version,
}: {
  store: AuspiceStore;
  panelsWidth: number;
  visibleHeight: number;
  entropy: boolean;
  runs: FigureRuns;
  version: string | undefined;
}) {
  const panelsToDisplay = useAuspiceSelector(store, selectPanelsToDisplay);
  const panelLayout = useAuspiceSelector(store, selectPanelLayout);
  const twoTrees = useAuspiceSelector(store, selectTwoTrees);

  const shownPanels = useMemo(
    () => (entropy ? panelsToDisplay : panelsToDisplay.filter((panel) => panel !== "entropy")),
    [entropy, panelsToDisplay],
  );

  const panelRef = useRef<HTMLDivElement>(null);
  const zoomButtonsOffset = useZoomButtonsOffset(panelRef);
  const { ref: actionsRef, width: actionsWidth } = useElementSize<HTMLDivElement>();

  const layout = { panelsWidth, visibleHeight, panelsToDisplay: shownPanels, panelLayout, twoTrees };
  const overlaid = panelSizes({ ...layout, headerStacked: false });

  const stacked =
    twoTrees || overlaid.tree.width - zoomButtonsOffset < actionsWidth + HEADER_ITEM_GAP + HEADER_MIN_SUMMARY_WIDTH;

  const sizes = stacked ? panelSizes({ ...layout, headerStacked: true }) : overlaid;

  if (panelsWidth <= 0 || visibleHeight <= 0) {
    return null;
  }

  return (
    <div className="flex flex-wrap content-start justify-center gap-4 px-5 pb-4">
      <TreeCard
        store={store}
        panelRef={panelRef}
        actionsRef={actionsRef}
        size={sizes.tree}
        twoTrees={twoTrees}
        stacked={stacked}
        headerRight={zoomButtonsOffset}
        runs={runs}
        version={version}
      />
      <MapCard sizes={sizes} />
      <EntropyCard sizes={sizes} />
    </div>
  );
}

function TreeCard({
  store,
  panelRef,
  actionsRef,
  size,
  twoTrees,
  stacked,
  headerRight,
  runs,
  version,
}: {
  store: AuspiceStore;
  panelRef: RefObject<HTMLDivElement | null>;
  actionsRef: RefCallback<HTMLDivElement | null>;
  size: DrawingSize;
  twoTrees: boolean;
  stacked: boolean;
  headerRight: number;
  runs: FigureRuns;
  version: string | undefined;
}) {
  const treeKey = useAuspiceSelector(store, selectSecondTreeName);

  const style = useMemo(
    () => cardStyle(twoTrees ? CARD_TITLE_HEIGHT_TWO_TREES : CARD_TITLE_HEIGHT, headerRight),
    [headerRight, twoTrees],
  );

  return (
    <div
      ref={panelRef}
      className="auspice-card relative flex-none"
      data-header={stacked ? "stacked" : "overlay"}
      // oxlint-disable-next-line react/forbid-dom-props -- the title row height and the zoom button position are runtime values that no Tailwind class can hold
      style={style}
    >
      <div
        className={
          stacked
            ? "relative flex h-[26px] items-center gap-4"
            : "absolute top-0 right-(--auspice-header-right) left-0 z-10 flex h-[26px] items-center gap-4"
        }
      >
        <div className="flex min-w-0 flex-1 items-center overflow-hidden whitespace-nowrap">
          <FiltersSummary />
        </div>
        <div ref={actionsRef} className="flex flex-none items-center gap-4 self-stretch">
          <TreeDownloads store={store} runs={runs} version={version} />
        </div>
      </div>
      {size.height > 0 && <Tree key={treeKey} width={size.width} height={size.height} />}
    </div>
  );
}

function MapCard({ sizes }: { sizes: PanelSizes }) {
  const webGl = useMemo(() => canDrawMaps(), []);
  const legend = useMemo(() => ({ swatches: !sizes.mapInGrid, demes: true }), [sizes.mapInGrid]);

  if (sizes.map === undefined || sizes.map.height <= 0) {
    return null;
  }

  if (!webGl) {
    return (
      <Alert className="basis-full">
        <AlertTitle>No map</AlertTitle>
        <AlertDescription>
          This browser cannot draw the map, because WebGL is turned off or not supported.
        </AlertDescription>
      </Alert>
    );
  }

  return (
    <div className={sizes.mapInGrid ? "auspice-card relative flex-none" : "auspice-card relative flex-none basis-full"}>
      <ErrorBoundary FallbackComponent={MapFailure}>
        <Suspense fallback={null}>
          <AuspiceMap
            width={sizes.map.width}
            height={sizes.map.height}
            justGotNewDatasetRenderNewMap={false}
            legend={legend}
          />
        </Suspense>
      </ErrorBoundary>
    </div>
  );
}

function MapFailure({ error }: FallbackProps) {
  return (
    <Alert variant="destructive">
      <AlertTitle>Auspice cannot draw the map</AlertTitle>
      <AlertDescription>{errorMessage(error)}</AlertDescription>
    </Alert>
  );
}

function EntropyCard({ sizes }: { sizes: PanelSizes }) {
  if (sizes.entropy === undefined || sizes.entropy.height <= 0) {
    return null;
  }

  return (
    <div className="auspice-card relative flex-none basis-full select-none">
      <Entropy width={sizes.entropy.width} height={sizes.entropy.height} />
    </div>
  );
}

function useZoomButtonsOffset(panelRef: RefObject<HTMLDivElement | null>): number {
  const [offset, setOffset] = useState(0);
  const [watcher] = useState(() => new ZoomButtonsWatcher(setOffset));

  useEffect(() => () => watcher.stop(), [watcher]);

  useLayoutEffect(() => {
    watcher.watch(panelRef.current);
  });

  return offset;
}

class ZoomButtonsWatcher {
  readonly #observer: ResizeObserver;
  #panel: HTMLElement | null = null;
  #buttons: Element | null = null;

  constructor(onOffset: (offset: number) => void) {
    this.#observer = new ResizeObserver(() => onOffset(this.#offset()));
  }

  watch(panel: HTMLElement | null): void {
    const buttons = panel?.querySelector(ZOOM_BUTTONS) ?? null;

    if (panel === this.#panel && buttons === this.#buttons) {
      return;
    }

    this.stop();
    this.#panel = panel;
    this.#buttons = buttons;

    if (panel !== null) {
      this.#observer.observe(panel);
    }

    if (buttons !== null) {
      this.#observer.observe(buttons);
    }
  }

  stop(): void {
    this.#observer.disconnect();
    this.#panel = null;
    this.#buttons = null;
  }

  #offset(): number {
    if (this.#panel === null || this.#buttons === null) {
      return 0;
    }

    return this.#panel.getBoundingClientRect().right - this.#buttons.getBoundingClientRect().left + ZOOM_BUTTONS_GAP;
  }
}

function cardStyle(titleHeight: number, headerRight: number): CSSProperties {
  return {
    "--auspice-card-title-height": `${titleHeight}px`,
    "--auspice-header-right": `${headerRight}px`,
  };
}

function selectPanelsToDisplay(state: AuspiceState): readonly string[] {
  return state.controls.panelsToDisplay;
}

function selectPanelLayout(state: AuspiceState): AuspicePanelLayout {
  return state.controls.panelLayout;
}

function selectSecondTreeName(state: AuspiceState): string {
  return state.controls.showTreeToo === false ? "" : state.controls.showTreeToo;
}

function selectTwoTrees(state: AuspiceState): boolean {
  return state.controls.showTreeToo !== false;
}
