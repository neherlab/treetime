import { stopAction, useAppNoticesStore } from "../noticeStore";

beforeEach(() => {
  useAppNoticesStore.setState({ stop: undefined, stops: 0, chunkFailed: false });
});

describe("back end stop notice", () => {
  test.each([
    { name: "crash in the web app", restarts: true, desktop: false, expected: "dismiss" },
    { name: "crash in the desktop app", restarts: true, desktop: true, expected: "dismiss" },
    { name: "failed start in the web app", restarts: false, desktop: false, expected: "reload-page" },
    { name: "failed start in the desktop app", restarts: false, desktop: true, expected: "restart-back-end" },
  ])("$name offers $expected", ({ restarts, desktop, expected }) => {
    expect(stopAction(restarts, desktop)).toBe(expected);
  });

  test("a crash shows the first line of its reason and counts one stop", () => {
    const { stopped } = useAppNoticesStore.getState();

    stopped({
      reason: "the back end stopped and restarts: panicked at src/run.rs:1:1\nstack backtrace:",
      restarts: true,
      requested: false,
    });

    const { stop, stops } = useAppNoticesStore.getState();
    expect({ stop, stops }).toStrictEqual({
      stop: { message: "the back end stopped and restarts: panicked at src/run.rs:1:1", restarts: true },
      stops: 1,
    });
  });

  test("a failed start shows its start error", () => {
    const { stopped } = useAppNoticesStore.getState();

    stopped({
      reason: "the back end could not start: When opening the runs folder '/r': denied",
      restarts: false,
      requested: false,
      error: { code: "internal_error", message: "When opening the runs folder '/r'", causes: ["denied"] },
    });

    expect(useAppNoticesStore.getState().stop).toStrictEqual({
      message: "When opening the runs folder '/r': denied",
      restarts: false,
    });
  });

  test("a requested restart shows no notice", () => {
    const { stopped } = useAppNoticesStore.getState();

    stopped({ reason: "the back end restarts to open the new runs folder", restarts: true, requested: true });

    const { stop, stops } = useAppNoticesStore.getState();
    expect({ stop, stops }).toStrictEqual({ stop: undefined, stops: 0 });
  });

  test("dismiss hides the notice and keeps the count of stops", () => {
    const { stopped, dismissStop } = useAppNoticesStore.getState();
    stopped({ reason: "crashed", restarts: true, requested: false });

    dismissStop();

    const { stop, stops } = useAppNoticesStore.getState();
    expect({ stop, stops }).toStrictEqual({ stop: undefined, stops: 1 });
  });
});
