import type { AnyAction } from "redux";
import type { ThunkDispatch } from "redux-thunk";

import type {
  AuspiceEntropyState,
  AuspiceMetadataState,
  AuspiceNode,
  AuspicePublication,
  AuspiceState,
  AuspiceTranslate,
  AuspiceTreeState,
} from "../../../state";

type AuspiceDispatch = ThunkDispatch<AuspiceState, undefined, AnyAction>;

export declare function SVG(
  dispatch: AuspiceDispatch,
  t: AuspiceTranslate,
  metadata: AuspiceMetadataState,
  nodes: readonly AuspiceNode[],
  visibility: readonly number[],
  filePrefix: string,
  panelsInDOM: readonly string[],
  panelLayout: string,
  publications: readonly AuspicePublication[],
  caption?: readonly string[],
): Promise<void>;

export declare function exportTree(options: {
  dispatch: AuspiceDispatch;
  filePrefix: string;
  tree: AuspiceTreeState;
  isNewick?: boolean;
  temporal: boolean;
  colorings?: AuspiceMetadataState["colorings"];
  colorBy?: string;
}): void;

export declare function strainTSV(
  dispatch: AuspiceDispatch,
  filePrefix: string,
  nodes: readonly AuspiceNode[],
  nodeVisibilities: readonly number[],
): void;

export declare function entropyTSV(dispatch: AuspiceDispatch, filePrefix: string, entropy: AuspiceEntropyState): void;

export declare function auspiceJSON(dispatch: AuspiceDispatch, state: AuspiceState, filePrefix: string): void;
