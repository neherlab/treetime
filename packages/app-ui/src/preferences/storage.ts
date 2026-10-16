import { zUiDraft, zUiSettings, type UiDraft, type UiSettings } from "@neherlab/app-contracts";
import { appSettings, appSettingsUi, type ApiClient } from "@neherlab/app-contracts/client";

export const PREFERENCES_STORAGE_KEY = "treetime-preferences";

export const DRAFT_STORAGE_KEY = "treetime-draft";

export function browserPreferencesStorage(
  shared: BrowserStorage,
  tab: BrowserStorage,
  reportInvalid: (message: string) => void = reportToConsole,
): PreferencesStorage {
  return new BrowserPreferencesStorage(shared, tab, reportInvalid);
}

export function apiPreferencesStorage(client: ApiClient): PreferencesStorage {
  return {
    load: async () => (await appSettings({ client, throwOnError: true })).data.ui ?? {},
    save: async (preferences) => {
      await appSettingsUi({ client, body: preferences, throwOnError: true });
    },
  };
}

export interface PreferencesStorage {
  load(): Promise<UiSettings>;
  save(preferences: UiSettings): Promise<void>;
}

export type BrowserStorage = Pick<Storage, "getItem" | "setItem">;

export function reportToConsole(message: string): void {
  console.error(`[TreeTime] ${message}`);
}

class BrowserPreferencesStorage implements PreferencesStorage {
  readonly #shared: BrowserStorage;
  readonly #tab: BrowserStorage;
  readonly #reportInvalid: (message: string) => void;
  #sharedText: string | null = null;

  constructor(shared: BrowserStorage, tab: BrowserStorage, reportInvalid: (message: string) => void) {
    this.#shared = shared;
    this.#tab = tab;
    this.#reportInvalid = reportInvalid;
  }

  load(): Promise<UiSettings> {
    const loaded = sharedPreferences(
      this.#parse(this.#shared.getItem(PREFERENCES_STORAGE_KEY), assertUiSettings, "preferences") ?? {},
    );

    this.#sharedText = JSON.stringify(loaded);

    const draft = this.#parse(this.#tab.getItem(DRAFT_STORAGE_KEY), assertUiDraft, "unfinished form");

    if (draft !== undefined) {
      loaded.draft = draft;
    }

    return Promise.resolve(loaded);
  }

  save({ draft, ...preferences }: UiSettings): Promise<void> {
    const text = JSON.stringify(sharedPreferences(preferences));

    if (text !== this.#sharedText) {
      this.#shared.setItem(PREFERENCES_STORAGE_KEY, text);
      this.#sharedText = text;
    }

    if (draft !== undefined) {
      this.#tab.setItem(DRAFT_STORAGE_KEY, JSON.stringify(draft));
    }

    return Promise.resolve();
  }

  #parse<T>(text: string | null, assert: (value: unknown) => asserts value is T, label: string): T | undefined {
    if (text === null) {
      return undefined;
    }

    try {
      const value: unknown = JSON.parse(text);

      assert(value);

      return value;
    } catch (error: unknown) {
      this.#reportInvalid(`ignored the ${label} saved in this browser, which cannot be read: ${String(error)}`);

      return undefined;
    }
  }
}

function sharedPreferences({ sidebar_width: sidebarWidth }: UiSettings): UiSettings {
  return sidebarWidth === undefined ? {} : { sidebar_width: sidebarWidth };
}

function assertUiSettings(value: unknown): asserts value is UiSettings {
  zUiSettings.parse(value);
}

function assertUiDraft(value: unknown): asserts value is UiDraft {
  zUiDraft.parse(value);
}
