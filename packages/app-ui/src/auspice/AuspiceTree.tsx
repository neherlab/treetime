import ChooseBranchLabelling from "auspice/src/components/controls/choose-branch-labelling";
import ChooseLayout from "auspice/src/components/controls/choose-layout";
import ChooseMetric from "auspice/src/components/controls/choose-metric";
import ChooseTipLabel from "auspice/src/components/controls/choose-tip-label";
import ColorBy, { ColorByInfo } from "auspice/src/components/controls/color-by";
import { ControlHeader } from "auspice/src/components/controls/controlHeader";
import FilterData, { FilterInfo } from "auspice/src/components/controls/filter";
import { TreeInfo } from "auspice/src/components/controls/miscInfoText";
import { ControlsContainer } from "auspice/src/components/controls/styles";
import { ToggleFocus } from "auspice/src/components/controls/toggle-focus";
import { DownloadButtons } from "auspice/src/components/download/downloadButtons";
import { publications } from "auspice/src/components/download/downloadModal";
import Entropy from "auspice/src/components/entropy";
import FiltersSummary from "auspice/src/components/info/filtersSummary";
import Tree from "auspice/src/components/tree";
import { calcUsableWidth } from "auspice/src/util/computeResponsive";
import { useCallback, useState } from "react";
import { I18nextProvider } from "react-i18next";
import { Provider } from "react-redux";
import { ThemeProvider } from "styled-components";

import { useElementWidth } from "../hooks/useElementWidth";
import { Button } from "../ui";
import { AUSPICE_I18N } from "./i18n";
import type { AuspiceState } from "./state";
import type { AuspiceStore } from "./store";
import { useAuspiceSelector } from "./store-hooks";

const SIDEBAR_THEME = {
  background: "var(--color-surface-2)",
  color: "var(--color-ink)",
  "font-family": "Lato, Helvetica Neue, Helvetica, sans-serif",
  sidebarBoxShadow: "rgba(0, 0, 0, 0.15)",
  selectedColor: "var(--color-accent)",
  unselectedColor: "var(--color-ink-muted)",
  alternateBackground: "var(--color-line-strong)",
};

const RELEVANT_PUBLICATIONS = [publications.treetime];

const MIN_TREE_WIDTH = 320;

const TREE_ROW_PX = 12;

const MIN_TREE_HEIGHT = 480;

const MAX_TREE_HEIGHT = 1100;

const ENTROPY_HEIGHT = 300;

export function AuspiceTree({ store, tips }: { store: AuspiceStore; tips: number }) {
  const [downloadsOpen, setDownloadsOpen] = useState(false);
  const toggleDownloads = useCallback(() => setDownloadsOpen((open) => !open), []);
  const showEntropy = useAuspiceSelector(store, selectShowEntropy);

  return (
    <I18nextProvider i18n={AUSPICE_I18N}>
      <ThemeProvider theme={SIDEBAR_THEME}>
        <Provider store={store}>
          <div className="light-scope border-line grid min-w-0 grid-cols-1 overflow-hidden rounded-lg border lg:grid-cols-[260px_minmax(0,1fr)]">
            <aside aria-label="Tree controls" className="border-line border-b lg:border-r lg:border-b-0">
              <ControlsContainer>
                <ControlHeader title="Color By" tooltip={ColorByInfo} />
                <ColorBy />
                <ControlHeader title="Filter Data" tooltip={FilterInfo} />
                <FilterData measurementsOn={false} />
                <ControlHeader title="Tree" tooltip={TreeInfo} />
                <ChooseLayout />
                <ChooseMetric />
                <ToggleFocus />
                <ChooseBranchLabelling />
                <ChooseTipLabel />
              </ControlsContainer>
            </aside>
            <div className="min-w-0">
              <div className="flex items-start gap-2 px-3 pt-2">
                <div className="min-w-0 flex-1">
                  <FiltersSummary />
                </div>
                <Button
                  type="button"
                  variant="outline"
                  size="sm"
                  onClick={toggleDownloads}
                  aria-expanded={downloadsOpen}
                >
                  {downloadsOpen ? "Hide downloads" : "Download figure and data"}
                </Button>
              </div>
              {downloadsOpen && (
                <div className="mx-3 mt-2 rounded-md bg-[#30353f] px-4 py-3">
                  <DownloadButtons relevantPublications={RELEVANT_PUBLICATIONS} />
                </div>
              )}
              <SizedPanels tips={tips} showEntropy={showEntropy} />
            </div>
          </div>
        </Provider>
      </ThemeProvider>
    </I18nextProvider>
  );
}

function SizedPanels({ tips, showEntropy }: { tips: number; showEntropy: boolean }) {
  const [container, setContainer] = useState<HTMLDivElement | null>(null);
  const width = Math.floor(calcUsableWidth(useElementWidth(container), 1));
  const height = Math.min(MAX_TREE_HEIGHT, Math.max(MIN_TREE_HEIGHT, tips * TREE_ROW_PX));

  return (
    <div ref={setContainer} className="relative min-w-0 pb-2">
      {width >= MIN_TREE_WIDTH && (
        <>
          <Tree width={width} height={height} />
          {showEntropy && <Entropy width={width} height={ENTROPY_HEIGHT} />}
        </>
      )}
    </div>
  );
}

function selectShowEntropy(state: AuspiceState): boolean {
  return state.controls.panelsToDisplay.includes("entropy");
}
