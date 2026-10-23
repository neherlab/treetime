export declare function calcPanelDims(
  panels: readonly string[],
  narrativeIsDisplayed: boolean,
  availableWidth: number,
  availableHeight: number,
): AuspicePanelDims;

interface AuspicePanelDims {
  full: AuspicePanelSize;
  grid: AuspicePanelSize;
  chartEntropy: AuspicePanelSize;
}

interface AuspicePanelSize {
  width: number;
  height: number;
}
