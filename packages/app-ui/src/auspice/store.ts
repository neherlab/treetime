import { performanceFlags } from "auspice/src/middleware/performanceFlags";
import { keepScatterplotStateInSync } from "auspice/src/middleware/scatterplot";
import browserDimensions from "auspice/src/reducers/browserDimensions";
import controls from "auspice/src/reducers/controls";
import entropy from "auspice/src/reducers/entropy";
import frequencies from "auspice/src/reducers/frequencies";
import measurements from "auspice/src/reducers/measurements";
import metadata from "auspice/src/reducers/metadata";
import narrative from "auspice/src/reducers/narrative";
import notifications from "auspice/src/reducers/notifications";
import tree from "auspice/src/reducers/tree";
import treeToo from "auspice/src/reducers/tree/treeToo";
import { applyMiddleware, combineReducers, legacy_createStore, type AnyAction, type Store } from "redux";
import thunk, { type ThunkDispatch } from "redux-thunk";

import type { AuspiceState } from "./state";

export type AuspiceStore = Store<AuspiceState> & {
  dispatch: ThunkDispatch<AuspiceState, undefined, AnyAction>;
};

const GENERAL_STATE = {
  defaults: { language: "en" },
  language: "en",
  mobileDisplay: false,
  displayComponent: "main",
  pathname: "",
};

export function createAuspiceStore(): AuspiceStore {
  const reducer = combineReducers({
    metadata,
    tree,
    frequencies,
    controls,
    entropy,
    browserDimensions,
    notifications,
    narrative,
    treeToo,
    measurements,
    general: () => GENERAL_STATE,
    query: () => ({}),
  });

  return legacy_createStore(
    reducer,
    applyMiddleware<ThunkDispatch<AuspiceState, undefined, AnyAction>, AuspiceState>(
      thunk,
      keepScatterplotStateInSync,
      performanceFlags,
    ),
  );
}
