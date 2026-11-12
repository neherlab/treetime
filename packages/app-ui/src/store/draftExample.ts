import type { UiDraft, UiDraftSource } from "@neherlab/app-contracts";

import { folderName } from "../settings/inputs";

export function withDraftExample(draft: UiDraft, previous: string | undefined): UiDraft {
  const { example: _example, ...rest } = draft;
  const example = draftExample(Object.values(draft.sources), previous);

  return example === undefined ? rest : { ...rest, example };
}

export function draftExample(sources: readonly UiDraftSource[], previous: string | undefined): string | undefined {
  if (sources.length === 0) {
    return undefined;
  }

  if (sources.every((source) => source.origin.kind === "run")) {
    return previous;
  }

  const folders = new Set(sources.map(exampleFolder));
  const [folder, ...others] = folders;

  return others.length === 0 ? folder : undefined;
}

function exampleFolder(source: UiDraftSource): string | undefined {
  if (source.origin.kind === "example") {
    return folderName(source.origin.path);
  }

  return source.origin.kind === "dataset" ? folderName(source.label) : undefined;
}
