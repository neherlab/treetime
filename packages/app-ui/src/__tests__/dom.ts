import { GlobalRegistrator } from "@happy-dom/global-registrator";
import { configure } from "@testing-library/dom";
import { act } from "react";

/* HACK: the DOM stays registered for the rest of the test file's worker (`bun test --parallel`), because React's
   scheduler runs a task queued before the last unmount that reads `window`, after an `unregister()` would remove it. */
export function registerDom() {
  GlobalRegistrator.register();
  // oxlint-disable-next-line typescript/strict-void-return -- user-event returns the wrapper's result as the result of `dispatchEvent`, although `@testing-library/dom` types the wrapper as void
  configure({ eventWrapper: inAct });
}

export function inAct(callback: () => unknown): unknown {
  let result: unknown;
  Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

  try {
    act(() => {
      result = callback();
    });
  } finally {
    Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: false });
  }

  return result;
}
