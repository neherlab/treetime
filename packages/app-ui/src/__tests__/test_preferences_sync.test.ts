import type { UiSettings } from "@neherlab/app-contracts";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";

import type { PreferencesStorage } from "../preferences/storage";
import { applyPreferences, currentPreferences, PreferencesSaver, restorePreferences } from "../preferences/sync";
import { useDraftStore } from "../store/draft";
import { freshDraft } from "../store/draftSchema";
import { usePreferencesStore } from "../store/preferences";
import { FakeServer, json } from "./api_server";

describe("sync preferences", () => {
  beforeEach(() => {
    usePreferencesStore.setState({ preferences: {} });
    useDraftStore.setState({ draft: freshDraft("timetree"), epoch: 0 });
  });

  test("loaded preferences reach the stores and the theme", () => {
    const themes: string[] = [];
    const draft = { ...freshDraft("clock"), search: "rate", from_run_id: "r1" };

    applyPreferences({ theme: "dark", sidebar_width: 400, draft }, (theme) => {
      themes.push(theme);
    });

    expect([themes, currentPreferences(), useDraftStore.getState().epoch]).toStrictEqual([
      ["dark"],
      { theme: "dark", sidebar_width: 400, draft },
      1,
    ]);
  });

  test("empty preferences keep the defaults and leave the theme alone", () => {
    const themes: string[] = [];

    applyPreferences({}, (theme) => {
      themes.push(theme);
    });

    expect([themes, currentPreferences()]).toStrictEqual([[], { draft: freshDraft("timetree") }]);
  });
});

describe("restore preferences", () => {
  const draft = {
    ...freshDraft("clock"),
    config: { ...freshDraft("clock").config, tree: "/treetime/runs/r1/inputs/tree.nwk" },
    sources: { tree: { label: "tree.nwk", origin: "upload" as const, run_id: "r1" } },
    upload_run_id: "r1",
  };

  test("a saved form loses the files of a run that the back end no longer has", async () => {
    const server = new FakeServer({});

    await expect(restorePreferences(server.client(), { theme: "dark", draft }, failTest)).resolves.toStrictEqual({
      preferences: { theme: "dark", draft: freshDraft("clock") },
      dropped: ["tree.nwk"],
    });
  });

  test("preferences without a form ask the back end nothing", async () => {
    const server = new FakeServer({});

    await expect(restorePreferences(server.client(), { sidebar_width: 300 }, failTest)).resolves.toStrictEqual({
      preferences: { sidebar_width: 300 },
      dropped: [],
    });
    expect(server.keys()).toStrictEqual([]);
  });

  test("a failed lookup keeps the saved form and is reported", async () => {
    const failures: string[] = [];
    const server = new FakeServer({
      "GET /api/runs/r1": () => json({ code: "internal_error", message: "the server is down", causes: [] }, 500),
    });

    const restored = await restorePreferences(server.client(), { draft }, (message) => {
      failures.push(message);
    });

    expect([restored, failures]).toStrictEqual([
      { preferences: { draft }, dropped: [] },
      ["the runs of the saved form could not be checked: the server is down"],
    ]);
  });
});

describe("sync saver", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    usePreferencesStore.setState({ preferences: {} });
    useDraftStore.setState({ draft: freshDraft("timetree"), epoch: 0 });
  });

  afterEach(() => {
    vi.useRealTimers();
  });

  test("a burst of changes saves the last preferences once", async () => {
    const storage = recordingStorage();
    const saver = new PreferencesSaver(storage, fakeTarget(), failTest);

    usePreferencesStore.getState().setSidebarWidth(300);
    usePreferencesStore.getState().setSidebarWidth(310);
    useDraftStore.getState().update({ search: "clock" });
    await vi.advanceTimersByTimeAsync(1_000);
    saver.dispose();

    expect(storage.saved).toStrictEqual([
      { sidebar_width: 310, draft: { ...freshDraft("timetree"), search: "clock" } },
    ]);
  });

  test("leaving the page saves a pending change at once", async () => {
    const storage = recordingStorage();
    const target = fakeTarget();
    const saver = new PreferencesSaver(storage, target, failTest);

    usePreferencesStore.getState().setTheme("light");
    target.emit("pagehide");
    await vi.advanceTimersByTimeAsync(0);
    saver.dispose();

    expect(storage.saved.map((preferences) => preferences.theme)).toStrictEqual(["light"]);
  });

  test("a disposed saver stops listening and saves its pending change", async () => {
    const storage = recordingStorage();
    const target = fakeTarget();
    const saver = new PreferencesSaver(storage, target, failTest);

    usePreferencesStore.getState().setTheme("dark");
    saver.dispose();
    usePreferencesStore.getState().setTheme("light");
    await vi.advanceTimersByTimeAsync(1_000);

    expect([storage.saved.map((preferences) => preferences.theme), target.listeners()]).toStrictEqual([["dark"], 0]);
  });

  test("saves run one after another in the order of the changes", async () => {
    const firstSave = Promise.withResolvers<undefined>();
    const storage = recordingStorage(firstSave.promise);
    const saver = new PreferencesSaver(storage, fakeTarget(), failTest);

    usePreferencesStore.getState().setTheme("dark");
    await vi.advanceTimersByTimeAsync(1_000);
    usePreferencesStore.getState().setTheme("light");
    await vi.advanceTimersByTimeAsync(1_000);
    firstSave.resolve(undefined);
    await vi.advanceTimersByTimeAsync(0);
    saver.dispose();

    expect(storage.saved.map((preferences) => preferences.theme)).toStrictEqual(["dark", "light"]);
  });

  test("a failed save is reported", async () => {
    const failures: string[] = [];

    const storage: PreferencesStorage = {
      load: () => Promise.resolve({}),
      save: () => Promise.reject(new Error("the back end stopped")),
    };

    const saver = new PreferencesSaver(storage, fakeTarget(), (message) => {
      failures.push(message);
    });

    usePreferencesStore.getState().setTheme("dark");
    await vi.advanceTimersByTimeAsync(1_000);
    saver.dispose();

    expect(failures).toStrictEqual(["the back end stopped"]);
  });
});

function recordingStorage(
  firstSave: Promise<undefined> = Promise.resolve(undefined),
): PreferencesStorage & { saved: UiSettings[] } {
  const saved: UiSettings[] = [];
  let calls = 0;

  return {
    saved,
    load: () => Promise.resolve({}),
    save: async (preferences) => {
      calls += 1;

      if (calls === 1) {
        await firstSave;
      }

      saved.push(preferences);
    },
  };
}

function fakeTarget() {
  const handlers = new Map<string, Set<() => void>>();

  return {
    addEventListener(type: string, listener: () => void) {
      handlers.set(type, (handlers.get(type) ?? new Set()).add(listener));
    },
    removeEventListener(type: string, listener: () => void) {
      handlers.get(type)?.delete(listener);
    },
    emit(type: string) {
      handlers.get(type)?.forEach((listener) => {
        listener();
      });
    },
    listeners: () => [...handlers.values()].reduce((count, set) => count + set.size, 0),
  };
}

function failTest(message: string): never {
  throw new Error(message);
}
