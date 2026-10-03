import { zUiSettings, type UiSettings } from "@neherlab/app-contracts";
import { appSettings, appSettingsUi, type ApiClient } from "@neherlab/app-contracts/client";
import type * as z from "zod";

export const PREFERENCES_STORAGE_KEY = "treetime-preferences";

export function browserPreferencesStorage(
  storage: Pick<Storage, "getItem" | "setItem">,
  reportInvalid: (message: string) => void = reportToConsole,
): PreferencesStorage {
  return {
    load: () => Promise.resolve(storedPreferences(storage.getItem(PREFERENCES_STORAGE_KEY), reportInvalid)),
    save: (preferences) => {
      storage.setItem(PREFERENCES_STORAGE_KEY, JSON.stringify(preferences));

      return Promise.resolve();
    },
  };
}

export function apiPreferencesStorage(client: ApiClient): PreferencesStorage {
  return {
    load: async () => (await appSettings({ client, throwOnError: true })).data.ui ?? {},
    save: async (preferences) => {
      await appSettingsUi({ client, body: preferences, throwOnError: true });
    },
  };
}

export type LoadedPreferences = z.output<typeof zUiSettings>;

export interface PreferencesStorage {
  load(): Promise<LoadedPreferences>;
  save(preferences: UiSettings): Promise<void>;
}

function storedPreferences(text: string | null, reportInvalid: (message: string) => void): LoadedPreferences {
  if (text === null) {
    return {};
  }

  try {
    const value: unknown = JSON.parse(text);

    return zUiSettings.parse(value);
  } catch (error: unknown) {
    reportInvalid(`the preferences saved in this browser are not valid and are ignored: ${String(error)}`);

    return {};
  }
}

function reportToConsole(message: string): void {
  console.error(`[TreeTime] ${message}`);
}
