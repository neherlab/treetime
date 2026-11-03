import { numberOfGridPanels } from "auspice/src/actions/panelDisplay";
import { calcPanelDims } from "auspice/src/components/main/utils";
import { controlsHiddenWidth, controlsPadding, controlsWidth, twoColumnBreakpoint } from "auspice/src/util/globals";

import type { AuspicePanelLayout } from "./state";

export const SIDEBAR_WIDTH = controlsWidth + controlsPadding;

export const CARD_TITLE_HEIGHT = 26;

export const CARD_TITLE_HEIGHT_TWO_TREES = 46;

const VIEW_PANELS = ["tree", "map", "entropy"] as const;

export function isSidebarDocked(viewerWidth: number): boolean {
  return viewerWidth > controlsHiddenWidth;
}

export function startLayout(viewerWidth: number): StartLayout {
  return {
    sidebarOpen: viewerWidth > controlsHiddenWidth,
    panelLayout: viewerWidth > twoColumnBreakpoint ? "grid" : "full",
  };
}

export interface StartLayout {
  sidebarOpen: boolean;
  panelLayout: AuspicePanelLayout;
}

export function panelSizes({
  panelsWidth,
  visibleHeight,
  panelsToDisplay,
  panelLayout,
  twoTrees,
  headerStacked,
}: PanelSizesParams): PanelSizes {
  const { full, grid, chartEntropy } = calcPanelDims(panelsToDisplay, false, panelsWidth, visibleHeight);
  const inGrid = panelLayout === "grid";
  const showMap = panelsToDisplay.includes("map") && !twoTrees;
  const mapInGrid = inGrid && numberOfGridPanels(panelsToDisplay) % 2 === 0;
  const treeCard = inGrid && !twoTrees ? grid : full;
  const mapCard = mapInGrid ? grid : full;

  const treeTitle =
    (twoTrees ? CARD_TITLE_HEIGHT_TWO_TREES : CARD_TITLE_HEIGHT) + (headerStacked ? CARD_TITLE_HEIGHT : 0);

  return {
    tree: drawingSize(treeCard, treeTitle),
    map: showMap ? drawingSize(mapCard, CARD_TITLE_HEIGHT) : undefined,
    mapInGrid,
    entropy: panelsToDisplay.includes("entropy") ? drawingSize(chartEntropy, CARD_TITLE_HEIGHT) : undefined,
  };
}

export function renderedPanels({
  panelsToDisplay,
  sizes,
  webGl,
  failed,
}: {
  panelsToDisplay: readonly string[];
  sizes: PanelSizes;
  webGl: boolean;
  failed: ReadonlySet<ViewPanel>;
}): readonly string[] {
  const drawn: Readonly<Record<ViewPanel, boolean>> = {
    tree: sizes.tree.height > 0,
    map: webGl && (sizes.map?.height ?? 0) > 0,
    entropy: (sizes.entropy?.height ?? 0) > 0,
  };

  return panelsToDisplay.filter((panel) => isViewPanel(panel) && drawn[panel] && !failed.has(panel));
}

export type ViewPanel = (typeof VIEW_PANELS)[number];

export interface PanelSizesParams {
  panelsWidth: number;
  visibleHeight: number;
  panelsToDisplay: readonly string[];
  panelLayout: AuspicePanelLayout;
  twoTrees: boolean;
  headerStacked: boolean;
}

export interface PanelSizes {
  tree: DrawingSize;
  map: DrawingSize | undefined;
  mapInGrid: boolean;
  entropy: DrawingSize | undefined;
}

export interface DrawingSize {
  width: number;
  height: number;
}

function drawingSize(card: DrawingSize, titleHeight: number): DrawingSize {
  return { width: Math.max(0, Math.floor(card.width)), height: Math.max(0, Math.floor(card.height - titleHeight)) };
}

function isViewPanel(panel: string): panel is ViewPanel {
  return VIEW_PANELS.some((viewPanel) => viewPanel === panel);
}
