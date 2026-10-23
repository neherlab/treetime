import type { Dispatch } from "redux";

import type { AuspiceState } from "../../state";

export declare function createStateFromQueryOrJSONs(params: {
  json: unknown;
  query: Readonly<Record<string, string>>;
  dispatch: Dispatch;
}): AuspiceCleanState;

export declare function createTreeTooState(params: {
  json: unknown;
  oldState: AuspiceState;
  originalTreeUrl: string;
  secondTreeUrl: string;
  dispatch: Dispatch;
}): AuspiceTreeTooState;

interface AuspiceCleanState {
  metadata: unknown;
  tree: unknown;
  treeToo: unknown;
  controls: unknown;
  entropy: unknown;
  frequencies: unknown;
  narrative: unknown;
  measurements: unknown;
  query: unknown;
}

interface AuspiceTreeTooState {
  tree: unknown;
  treeToo: unknown;
  controls: unknown;
  metadata: unknown;
}
