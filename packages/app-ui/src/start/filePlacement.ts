import { type DetectedFile, type InputKind, type InputSlot, zInputKind } from "@neherlab/app-contracts";

import { baseName } from "../settings/inputs";

const INPUT_KINDS: ReadonlySet<string> = new Set(zInputKind.options);

export interface Placement {
  inputs: PlacedInput[];
  configs: number[];
  packages: number[];
  unplaced: UnplacedFile[];
}

interface PlacedInput {
  key: InputKind;
  list: boolean;
  files: number[];
}

interface UnplacedFile {
  index: number;
  reason: string;
}

export function placeFiles(
  detected: readonly DetectedFile[],
  slots: readonly InputSlot[],
  analysis: string,
): Placement {
  const indexed = detected.map((file, index) => ({ file, index }));
  const ofKind = (kind: DetectedFile["kind"]) => indexed.filter(({ file }) => file.kind === kind);
  const [config, ...surplusConfigs] = ofKind("configuration");

  const slotted = slots.flatMap((slot) => {
    const files = ofKind(slot.kind);
    const kept = slot.list ? files : files.slice(0, 1);
    const [first] = kept;

    return first === undefined
      ? []
      : [
          {
            input: { key: slot.kind, list: slot.list, files: kept.map(({ index }) => index) },
            surplus: files.slice(kept.length).map(({ index }) => ({
              index,
              reason: `The analysis reads one ${slot.kind}; ${baseName(first.file.path)} fills it`,
            })),
          },
        ];
  });

  const readable = new Set<DetectedFile["kind"]>(slots.map((slot) => slot.kind));

  const unread = indexed.flatMap(({ file, index }) => {
    if (file.kind === "unknown") {
      return [{ index, reason: file.reason ?? "No reader of TreeTime reads the file" }];
    }

    return isInputKind(file.kind) && !readable.has(file.kind)
      ? [{ index, reason: `The ${analysis} analysis reads no ${file.kind}` }]
      : [];
  });

  const surplus = surplusConfigs.map(({ index }) => ({
    index,
    reason: `Only one configuration loads at a time; ${baseName(config?.file.path ?? "")} loads`,
  }));

  return {
    inputs: slotted.map(({ input }) => input),
    configs: config === undefined ? [] : [config.index],
    packages: ofKind("run-package").map(({ index }) => index),
    unplaced: [...slotted.flatMap((entry) => entry.surplus), ...unread, ...surplus].toSorted(
      (left, right) => left.index - right.index,
    ),
  };
}

function isInputKind(kind: DetectedFile["kind"]): kind is InputKind {
  return INPUT_KINDS.has(kind);
}
