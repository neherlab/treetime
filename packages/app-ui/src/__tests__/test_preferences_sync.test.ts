import type { UiSettings } from "@neherlab/app-contracts";
import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";

import type { PreferencesStorage } from "../preferences/storage";
import { applyPreferences, currentPreferences, PreferencesSaver } from "../preferences/sync";
import { useDraftStore } from "../store/draft";
import { freshDraft } from "../store/draftSchema";
import { usePreferencesStore } from "../store/preferences";

describe("sync preferences", () => {
  beforeEach(() => {
    usePreferencesStore.setState({ theme: undefined, sidebarWidth: undefined });
    useDraftStore.setState({ ...freshDraft("timetree"), epoch: 0 });
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

    expect([themes, currentPreferences()]).toStrictEqual([
      [],
      { theme: null, sidebar_width: null, draft: freshDraft("timetree") },
    ]);
  });
});

describe("sync saver", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    usePreferencesStore.setState({ theme: undefined, sidebarWidth: undefined });
    useDraftStore.setState({ ...freshDraft("timetree"), epoch: 0 });
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
      { theme: null, sidebar_width: 310, draft: { ...freshDraft("timetree"), search: "clock" } },
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
