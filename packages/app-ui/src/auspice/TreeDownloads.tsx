import { errorMessage } from "@neherlab/app-contracts";
import { publications } from "auspice/src/components/download/downloadModal";
import { auspiceJSON, entropyTSV, exportTree, strainTSV, SVG } from "auspice/src/components/download/helperFunctions";
import { datasetSummary } from "auspice/src/components/info/datasetSummary";
import { getMapTilesSettings } from "auspice/src/util/globals";
import { useCallback, useMemo } from "react";
import ChevronDownIcon from "~icons/lucide/chevron-down";

import { Button } from "../ui/button";
import { DropdownMenu, DropdownMenuContent, DropdownMenuItem, DropdownMenuTrigger } from "../ui/dropdown-menu";
import { useToastManager } from "../ui/toast";
import { figureCaption } from "./caption";
import { downloadPrefix, type FigureRuns } from "./downloads";
import { translateAuspice } from "./i18n";
import type { AuspiceSharing, AuspiceState } from "./state";
import type { AuspiceStore } from "./store";
import { useAuspiceSelector } from "./store-hooks";

const RELEVANT_PUBLICATIONS = [publications.treetime];

const TREE_DOWNLOADS: readonly TreeDownload[] = [
  {
    name: () => "Figure (SVG)",
    description: "The panels as shown, with the run named below",
    available: ({ sharing }) => sharing.screenshot === true,
    save: saveFigure,
  },
  {
    name: ({ temporal }) => `${treeKind(temporal)} (Newick)`,
    description: "The tree in view, with branch lengths of the shown metric",
    available: ({ sharing }) => sharing.trees === true,
    save: ({ store, prefix }) => {
      const state = store.getState();

      exportTree({
        isNewick: true,
        dispatch: store.dispatch,
        filePrefix: prefix,
        tree: state.tree,
        temporal: isTemporal(state),
      });
    },
  },
  {
    name: ({ temporal }) => `${treeKind(temporal)} (Nexus)`,
    description: "The tree in view, with the colorings as annotations",
    available: ({ sharing }) => sharing.trees === true,
    save: ({ store, prefix }) => {
      const state = store.getState();

      exportTree({
        dispatch: store.dispatch,
        filePrefix: prefix,
        tree: state.tree,
        colorings: state.metadata.colorings,
        colorBy: state.controls.colorBy,
        temporal: isTemporal(state),
      });
    },
  },
  {
    name: () => "Metadata (TSV)",
    description: "The attributes of the shown tips",
    available: ({ sharing }) => sharing.metadata_tsv === true,
    save: ({ store, prefix }) => {
      const state = store.getState();

      strainTSV(store.dispatch, prefix, state.tree.nodes ?? [], state.tree.visibility ?? []);
    },
  },
  {
    name: () => "Diversity (TSV)",
    description: "The values of the entropy chart",
    available: ({ sharing, entropyLoaded }) => sharing.entropy === true && entropyLoaded,
    save: ({ store, prefix }) => entropyTSV(store.dispatch, prefix, store.getState().entropy),
  },
  {
    name: () => "Auspice JSON",
    description: "The dataset with the current view as its default",
    available: ({ sharing }) => sharing.dataset_json === true,
    save: ({ store, prefix }) => auspiceJSON(store.dispatch, store.getState(), prefix),
  },
];

export function TreeDownloads({
  store,
  runs,
  version,
}: {
  store: AuspiceStore;
  runs: FigureRuns;
  version: string | undefined;
}) {
  const toasts = useToastManager();
  const sharing = useAuspiceSelector(store, selectSharing);
  const entropyLoaded = useAuspiceSelector(store, selectEntropyLoaded);
  const temporal = useAuspiceSelector(store, isTemporal);
  const flags = useMemo(() => ({ sharing, entropyLoaded, temporal }), [entropyLoaded, sharing, temporal]);
  const prefix = downloadPrefix(runs);

  const context = useMemo<DownloadContext>(
    () => ({
      store,
      prefix,
      runs,
      version,
      fail: (title, description) => {
        toasts.add({ title, description, type: "error", timeout: 0 });
      },
    }),
    [prefix, runs, store, toasts, version],
  );

  const downloads = TREE_DOWNLOADS.filter((download) => download.available(flags));

  return (
    <DropdownMenu>
      <DropdownMenuTrigger render={<Button type="button" variant="outline" size="xs" />}>
        Download
        <ChevronDownIcon aria-hidden data-icon="inline-end" />
      </DropdownMenuTrigger>
      <DropdownMenuContent align="end" className="w-72">
        {downloads.map((download) => (
          <DownloadItem key={download.description} download={download} flags={flags} context={context} />
        ))}
      </DropdownMenuContent>
    </DropdownMenu>
  );
}

interface TreeDownload {
  name: (flags: DownloadFlags) => string;
  description: string;
  available: (flags: DownloadFlags) => boolean;
  save: (context: DownloadContext) => void;
}

interface DownloadContext {
  store: AuspiceStore;
  prefix: string;
  runs: FigureRuns;
  version: string | undefined;
  fail: (title: string, description: string) => void;
}

interface DownloadFlags {
  sharing: AuspiceSharing;
  entropyLoaded: boolean;
  temporal: boolean;
}

function DownloadItem({
  download,
  flags,
  context,
}: {
  download: TreeDownload;
  flags: DownloadFlags;
  context: DownloadContext;
}) {
  const save = useCallback(() => download.save(context), [context, download]);

  return (
    <DropdownMenuItem onClick={save} className="flex-col items-start gap-0.5">
      <span>{download.name(flags)}</span>
      <span className="text-muted-foreground text-xs">{download.description}</span>
    </DropdownMenuItem>
  );
}

function saveFigure(context: DownloadContext): void {
  void writeFigure(context);
}

async function writeFigure({ store, prefix, runs, version, fail }: DownloadContext): Promise<void> {
  const { metadata, tree, controls } = store.getState();
  const nodes = tree.nodes ?? [];
  const visibility = tree.visibility ?? [];

  const caption = figureCaption({
    runs,
    version,
    updated: metadata.updated,
    summary: datasetSummary({
      nodes,
      visibility,
      mainTreeNumTips: metadata.mainTreeNumTips,
      branchLengthsToDisplay: controls.branchLengthsToDisplay,
      t: translateAuspice,
    }),
    mapAttribution: controls.panelsToDisplay.includes("map") ? getMapTilesSettings().attribution : undefined,
    publications: RELEVANT_PUBLICATIONS,
  });

  try {
    await SVG(
      store.dispatch,
      translateAuspice,
      metadata,
      nodes,
      visibility,
      prefix,
      controls.panelsToDisplay,
      controls.panelLayout,
      RELEVANT_PUBLICATIONS,
      caption,
    );
  } catch (error) {
    fail("The figure cannot be saved", errorMessage(error));
  }
}

function selectSharing(state: AuspiceState): AuspiceSharing {
  return state.metadata.sharing;
}

function selectEntropyLoaded(state: AuspiceState): boolean {
  return state.entropy.loaded;
}

function treeKind(temporal: boolean): string {
  return temporal ? "Time tree" : "Tree";
}

function isTemporal(state: Pick<AuspiceState, "controls">): boolean {
  return state.controls.distanceMeasure === "num_date";
}
