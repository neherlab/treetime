import type { UiDraft, UiDraftSource } from "@neherlab/app-contracts";

import { commandSettings } from "../settings/catalog";
import { folderName } from "../settings/inputs";

export function datasetLabel(draft: UiDraft): string | undefined {
  const order = commandSettings(draft.command).inputs.map((slot) => slot.kind);

  const sources = Object.entries(draft.sources)
    .toSorted(([left], [right]) => rank(order, left) - rank(order, right))
    .map(([, source]) => source);

  if (sources.length === 0) {
    return undefined;
  }

  const folders = new Set(sources.map(exampleFolder));
  const [folder] = folders;

  return folders.size === 1 && folder !== undefined ? folder : sources.map((source) => source.label).join(", ");
}

export function exampleDatasetName(draft: UiDraft): string | undefined {
  const folders = new Set(Object.values(draft.sources).map(exampleFolder));
  const [folder, ...others] = folders;

  return others.length === 0 ? folder : undefined;
}

function exampleFolder(source: UiDraftSource): string | undefined {
  if (source.origin.kind === "example") {
    return folderName(source.origin.path);
  }

  return source.origin.kind === "dataset" ? folderName(source.label) : undefined;
}

function rank(order: readonly string[], key: string): number {
  const index = order.indexOf(key);

  return index < 0 ? order.length : index;
}
