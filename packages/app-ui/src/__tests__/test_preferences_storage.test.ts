import type { UiSettings } from "@neherlab/app-contracts";

import {
  browserPreferencesStorage,
  type BrowserStorage,
  DRAFT_STORAGE_KEY,
  PREFERENCES_STORAGE_KEY,
  type PreferencesStorage,
} from "../preferences/storage";
import { freshDraft } from "../store/draftSchema";

describe("preferences storage browser", () => {
  const draft = { ...freshDraft("clock"), search: "rate" };

  test("nothing stored loads no preferences", async () => {
    await expect(new Browser().tab().load()).resolves.toStrictEqual({});
  });

  test("the sidebar width and the form load back after a reload", async () => {
    const browser = new Browser();
    const tab = browser.tab();

    await tab.save({ sidebar_width: 420, draft });

    await expect(tab.reload().load()).resolves.toStrictEqual({ sidebar_width: 420, draft });
  });

  test("the theme is not stored, because next-themes stores it and keeps the tabs in step", async () => {
    const browser = new Browser();

    await browser.tab().save({ theme: "dark", sidebar_width: 420 });

    expect(browser.shared.values).toStrictEqual(new Map([[PREFERENCES_STORAGE_KEY, '{"sidebar_width":420}']]));
  });

  test("a new tab starts with an empty form", async () => {
    const browser = new Browser();

    await browser.tab().save({ sidebar_width: 420, draft });

    await expect(browser.tab().load()).resolves.toStrictEqual({ sidebar_width: 420 });
  });

  test("a reload keeps the form of its own tab after another tab saved a form", async () => {
    const browser = new Browser();
    const first = browser.tab();
    const second = browser.tab();

    await first.save({ draft });
    await second.save({ draft: freshDraft("timetree") });

    await expect(first.reload().load()).resolves.toStrictEqual({ draft });
  });

  test("a form change in another tab keeps the sidebar width that this tab chose", async () => {
    const browser = new Browser();
    const first = browser.tab();
    const second = browser.tab();

    const loaded = await second.load();

    await first.load();
    await first.save({ sidebar_width: 500 });
    await second.save({ ...loaded, draft });

    await expect(first.reload().load()).resolves.toStrictEqual({ sidebar_width: 500 });
  });

  test("a sidebar width chosen in another tab replaces the stored one", async () => {
    const browser = new Browser();
    const first = browser.tab();
    const second = browser.tab();

    await first.load();
    await second.load();
    await first.save({ sidebar_width: 500 });
    await second.save({ sidebar_width: 300 });

    await expect(first.reload().load()).resolves.toStrictEqual({ sidebar_width: 300 });
  });

  test("preferences of an unknown shape are ignored and reported, and the form still loads", async () => {
    const reports: string[] = [];
    const browser = new Browser({ [PREFERENCES_STORAGE_KEY]: JSON.stringify({ theme: "blue" }) });

    const tab = browser.tab(
      (message) => {
        reports.push(message);
      },
      { [DRAFT_STORAGE_KEY]: JSON.stringify(draft) },
    );

    await expect(tab.load()).resolves.toStrictEqual({ draft });
    expect(reports).toHaveLength(1);
  });

  test("a form that is not JSON is ignored and reported, and the preferences still load", async () => {
    const reports: string[] = [];
    const browser = new Browser({ [PREFERENCES_STORAGE_KEY]: '{"sidebar_width":420}' });

    const tab = browser.tab(
      (message) => {
        reports.push(message);
      },
      { [DRAFT_STORAGE_KEY]: "{" },
    );

    await expect(tab.load()).resolves.toStrictEqual({ sidebar_width: 420 });
    expect(reports).toHaveLength(1);
  });
});

class Browser {
  readonly shared: MemoryStorage;

  constructor(shared: Record<string, string> = {}) {
    this.shared = new MemoryStorage(shared);
  }

  tab(reportInvalid: (message: string) => void = failTest, session: Record<string, string> = {}): Tab {
    return new Tab(this.shared, new MemoryStorage(session), reportInvalid);
  }
}

class Tab {
  readonly #shared: MemoryStorage;
  readonly #session: MemoryStorage;
  readonly #reportInvalid: (message: string) => void;
  readonly #storage: PreferencesStorage;

  constructor(shared: MemoryStorage, session: MemoryStorage, reportInvalid: (message: string) => void) {
    this.#shared = shared;
    this.#session = session;
    this.#reportInvalid = reportInvalid;
    this.#storage = browserPreferencesStorage(shared, session, reportInvalid);
  }

  async load(): Promise<UiSettings> {
    return this.#storage.load();
  }

  async save(preferences: UiSettings): Promise<void> {
    return this.#storage.save(preferences);
  }

  reload(): Tab {
    return new Tab(this.#shared, this.#session, this.#reportInvalid);
  }
}

class MemoryStorage implements BrowserStorage {
  readonly values: Map<string, string>;

  constructor(initial: Record<string, string>) {
    this.values = new Map(Object.entries(initial));
  }

  getItem(key: string): string | null {
    return this.values.get(key) ?? null;
  }

  setItem(key: string, value: string): void {
    this.values.set(key, value);
  }
}

function failTest(message: string): never {
  throw new Error(message);
}
