import ChooseBranchLabelling from "auspice/src/components/controls/choose-branch-labelling";
import ChooseExplodeAttr from "auspice/src/components/controls/choose-explode-attr";
import ChooseLayout from "auspice/src/components/controls/choose-layout";
import ChooseMetric from "auspice/src/components/controls/choose-metric";
import ChooseTipLabel from "auspice/src/components/controls/choose-tip-label";
import ColorBy, { ColorByInfo } from "auspice/src/components/controls/color-by";
import { ControlHeader } from "auspice/src/components/controls/controlHeader";
import DateRangeInputs, { DateRangeInfo } from "auspice/src/components/controls/date-range-inputs";
import FilterData, { FilterInfo } from "auspice/src/components/controls/filter";
import { EntropyInfo, ExplodeTreeInfo, PanelLayoutInfo, TreeInfo } from "auspice/src/components/controls/miscInfoText";
import PanelLayout from "auspice/src/components/controls/panel-layout";
import { PanelSection } from "auspice/src/components/controls/panelSection";
import { ControlsContainer } from "auspice/src/components/controls/styles";
import { ToggleFocus } from "auspice/src/components/controls/toggle-focus";
import type { ComponentProps } from "react";
import ChevronLeftIcon from "~icons/lucide/chevron-left";
import ChevronRightIcon from "~icons/lucide/chevron-right";

import { Button } from "../ui/button";
import { Tooltip, TooltipContent, TooltipTrigger } from "../ui/tooltip";
import type { AuspiceState } from "./state";
import type { AuspiceStore } from "./store";
import { useAuspiceSelector } from "./store-hooks";

const HIDE_CONTROLS = "Hide Auspice controls";

const SHOW_CONTROLS = "Show Auspice controls";

export function SidebarControls({ store, entropy }: { store: AuspiceStore; entropy: boolean }) {
  const dated = useAuspiceSelector(store, selectDated);
  const entropyAvailable = useAuspiceSelector(store, selectEntropyAvailable);
  const canTogglePanelLayout = useAuspiceSelector(store, selectCanTogglePanelLayout);

  return (
    <ControlsContainer>
      {dated && (
        <>
          <ControlHeader title="Date Range" tooltip={DateRangeInfo} />
          <DateRangeInputs />
        </>
      )}
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
      <ChooseExplodeAttr tooltip={ExplodeTreeInfo} />
      {entropy && entropyAvailable && <PanelSection panel="entropy" title="Entropy" tooltip={EntropyInfo} />}
      {canTogglePanelLayout && (
        <>
          <ControlHeader title="Display" tooltip={PanelLayoutInfo} />
          <PanelLayout />
        </>
      )}
    </ControlsContainer>
  );
}

export function HideControlsButton(props: ComponentProps<typeof Button>) {
  return (
    <Tooltip>
      <TooltipTrigger
        render={
          <Button
            type="button"
            variant="ghost"
            size="icon-xs"
            aria-label={HIDE_CONTROLS}
            className="absolute top-1 right-1 z-10 text-[#444]"
            {...props}
          />
        }
      >
        <ChevronLeftIcon aria-hidden className="size-3" />
      </TooltipTrigger>
      <TooltipContent side="right">{HIDE_CONTROLS}</TooltipContent>
    </Tooltip>
  );
}

export function ShowControlsTab(props: ComponentProps<"button">) {
  return (
    <Tooltip>
      <TooltipTrigger
        render={
          <button
            type="button"
            aria-label={SHOW_CONTROLS}
            className="bg-muted focus-visible:ring-ring/50 absolute top-1 left-0 z-10 flex h-11 w-3.5 items-center justify-center rounded-r-md text-[#444] shadow-md outline-none focus-visible:ring-3"
            {...props}
          />
        }
      >
        <ChevronRightIcon aria-hidden className="size-3" />
      </TooltipTrigger>
      <TooltipContent side="right">{SHOW_CONTROLS}</TooltipContent>
    </Tooltip>
  );
}

function selectDated(state: AuspiceState): boolean {
  return state.controls.branchLengthsToDisplay !== "divOnly";
}

function selectEntropyAvailable(state: AuspiceState): boolean {
  return state.controls.panelsAvailable.includes("entropy");
}

function selectCanTogglePanelLayout(state: AuspiceState): boolean {
  return state.controls.canTogglePanelLayout;
}
