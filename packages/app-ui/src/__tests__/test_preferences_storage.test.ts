import { describe, expect, test } from "vitest";

import { browserPreferencesStorage, PREFERENCES_STORAGE_KEY } from "../preferences/storage";

describe("preferences storage browser", () => {
  test("nothing stored loads no preferences", async () => {
    await expect(browserPreferencesStorage(memoryStorage()).load()).resolves.toStrictEqual({});
  });

  test("saved preferences load back", async () => {
    const storage = browserPreferencesStorage(memoryStorage());
    const preferences = { theme: "dark", sidebar_width: 420, draft: null } as const;

    await storage.save(preferences);

    await expect(storage.load()).resolves.toStrictEqual(preferences);
  });

  test("preferences of an unknown shape are ignored and reported", async () => {
    const reports: string[] = [];
    const stored = memoryStorage({ [PREFERENCES_STORAGE_KEY]: JSON.stringify({ theme: "blue" }) });

    await expect(
      browserPreferencesStorage(stored, (message) => {
        reports.push(message);
      }).load(),
    ).resolves.toStrictEqual({});
    expect(reports).toHaveLength(1);
  });

  test("text that is not JSON is ignored and reported", async () => {
    const reports: string[] = [];
    const stored = memoryStorage({ [PREFERENCES_STORAGE_KEY]: "{" });

    await expect(
      browserPreferencesStorage(stored, (message) => {
        reports.push(message);
      }).load(),
    ).resolves.toStrictEqual({});
    expect(reports).toHaveLength(1);
  });
});

function memoryStorage(initial: Record<string, string> = {}): Pick<Storage, "getItem" | "setItem"> {
  const values = new Map(Object.entries(initial));

  return {
    getItem: (key) => values.get(key) ?? null,
    setItem: (key, value) => {
      values.set(key, value);
    },
  };
}
