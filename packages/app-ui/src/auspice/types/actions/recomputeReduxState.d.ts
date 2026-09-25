import type { Dispatch } from "redux";

export declare function createStateFromQueryOrJSONs(params: {
  json: unknown;
  query: Readonly<Record<string, string>>;
  dispatch: Dispatch;
}): AuspiceCleanState;

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
