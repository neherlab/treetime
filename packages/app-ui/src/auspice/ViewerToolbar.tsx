import { Button } from "@neherlab/design/ui/button";
import {
  Combobox,
  ComboboxEmpty,
  ComboboxInput,
  ComboboxItem,
  ComboboxList,
  ComboboxPopup,
} from "@neherlab/design/ui/combobox";
import { NativeSelect, NativeSelectOption } from "@neherlab/design/ui/native-select";
import { changeLayout } from "auspice/src/actions/layout";
import { useCallback, useMemo, useState, type Ref } from "react";
import Maximize from "~icons/lucide/maximize";
import SlidersHorizontal from "~icons/lucide/sliders-horizontal";

import type { AuspiceState } from "./state";
import type { AuspiceStore } from "./store";
import { colorTreeBy, focusNode, showWholeTree, useAuspiceSelector } from "./store-hooks";
import { colorOptions, layoutOptions, tipMatches, tipNames } from "./toolbar";

export function ViewerToolbar({
  store,
  sidebarId,
  sidebarOpen,
  onSidebarOpen,
  optionsRef,
}: {
  store: AuspiceStore;
  sidebarId: string;
  sidebarOpen: boolean;
  onSidebarOpen: (open: boolean) => void;
  optionsRef: Ref<HTMLButtonElement>;
}) {
  const colorings = useAuspiceSelector(store, selectColorings);
  const colorBy = useAuspiceSelector(store, selectColorBy);
  const layout = useAuspiceSelector(store, selectLayout);
  const distanceMeasures = useAuspiceSelector(store, selectDistanceMeasures);
  const nodes = useAuspiceSelector(store, selectNodes);
  const zoomed = useAuspiceSelector(store, selectZoomed);
  const [query, setQuery] = useState("");

  const colors = useMemo(() => colorOptions(colorings, colorBy), [colorBy, colorings]);
  const layouts = useMemo(() => layoutOptions(distanceMeasures), [distanceMeasures]);
  const tips = useMemo(() => tipNames(nodes), [nodes]);
  const tipSet = useMemo(() => new Set(tips), [tips]);
  const matches = useMemo(() => tipMatches(tips, query), [query, tips]);

  const onColorBy = useCallback(
    (event: React.ChangeEvent<HTMLSelectElement>) => colorTreeBy(store, event.target.value),
    [store],
  );

  const onLayout = useCallback(
    (event: React.ChangeEvent<HTMLSelectElement>) => store.dispatch(changeLayout({ layout: event.target.value })),
    [store],
  );

  const onQuery = useCallback(
    (value: string) => {
      setQuery(value);

      if (tipSet.has(value)) {
        focusNode(store, value);
      }
    },
    [store, tipSet],
  );

  const onZoomOut = useCallback(() => showWholeTree(store), [store]);
  const onOptions = useCallback(() => onSidebarOpen(!sidebarOpen), [onSidebarOpen, sidebarOpen]);

  return (
    <fieldset aria-label="Tree view" className="flex min-w-0 flex-wrap items-center gap-1.5">
      <NativeSelect size="sm" aria-label="Color by" value={colorBy} onChange={onColorBy} className="w-44">
        {colors.map((option) => (
          <NativeSelectOption key={option.value} value={option.value}>
            {option.label}
          </NativeSelectOption>
        ))}
      </NativeSelect>
      <NativeSelect size="sm" aria-label="Tree layout" value={layout} onChange={onLayout} className="w-36">
        {layouts.map((option) => (
          <NativeSelectOption key={option.value} value={option.value}>
            {option.label}
          </NativeSelectOption>
        ))}
      </NativeSelect>
      <Combobox items={tips} filteredItems={matches} value={query} onValueChange={onQuery} autoHighlight>
        <ComboboxInput className="w-56" placeholder="Find a sample" aria-label="Find a sample in the tree" />
        <ComboboxPopup>
          <ComboboxEmpty>No sample matches</ComboboxEmpty>
          <ComboboxList>
            {(tip: string) => (
              <ComboboxItem key={tip} value={tip}>
                {tip}
              </ComboboxItem>
            )}
          </ComboboxList>
        </ComboboxPopup>
      </Combobox>
      <Button type="button" variant="ghost" size="sm" onClick={onZoomOut} disabled={!zoomed}>
        <Maximize aria-hidden />
        Zoom to root
      </Button>
      <Button
        ref={optionsRef}
        type="button"
        variant={sidebarOpen ? "secondary" : "outline"}
        size="sm"
        className="ml-auto"
        aria-expanded={sidebarOpen}
        aria-controls={sidebarId}
        onClick={onOptions}
      >
        <SlidersHorizontal aria-hidden />
        View options
      </Button>
    </fieldset>
  );
}

function selectColorings(state: AuspiceState): AuspiceState["metadata"]["colorings"] {
  return state.metadata.colorings;
}

function selectColorBy(state: AuspiceState): string {
  return state.controls.colorBy;
}

function selectLayout(state: AuspiceState): string {
  return state.controls.layout;
}

function selectDistanceMeasures(state: AuspiceState): readonly string[] {
  return state.controls.availableDistanceMeasures;
}

function selectNodes(state: AuspiceState): AuspiceState["tree"]["nodes"] {
  return state.tree.nodes;
}

function selectZoomed(state: AuspiceState): boolean {
  return state.tree.idxOfInViewRootNode !== 0;
}
