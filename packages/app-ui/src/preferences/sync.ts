import { errorMessage, type UiDraft, type UiSettings } from "@neherlab/app-contracts";
import { Debouncer } from "@tanstack/react-pacer";

import { useDraftStore } from "../store/draft";
import { storedDraft } from "../store/draftSchema";
import { usePreferencesStore } from "../store/preferences";
import type { LoadedPreferences, PreferencesStorage } from "./storage";

const SAVE_DELAY_MS = 300;

export function applyPreferences(preferences: LoadedPreferences, setTheme: (theme: string) => void): void {
  usePreferencesStore.setState({
    theme: preferences.theme ?? undefined,
    sidebarWidth: preferences.sidebar_width ?? undefined,
  });

  if (preferences.draft !== null && preferences.draft !== undefined) {
    useDraftStore.getState().load(storedDraft(preferences.draft));
  }

  if (preferences.theme !== null && preferences.theme !== undefined) {
    setTheme(preferences.theme);
  }
}

export function currentPreferences(): UiSettings {
  const { theme, sidebarWidth } = usePreferencesStore.getState();

  const { command, config, sources, from_run_id, upload_run_id, view, search, changed_only, code_format } =
    useDraftStore.getState();

  const draft: UiDraft = {
    command,
    config,
    sources,
    from_run_id,
    upload_run_id,
    view,
    search,
    changed_only,
    code_format,
  };

  return { theme: theme ?? null, sidebar_width: sidebarWidth ?? null, draft };
}

export class PreferencesSaver {
  readonly #storage: PreferencesStorage;
  readonly #target: PageEvents;
  readonly #onError: (message: string) => void;
  readonly #debouncer: Debouncer<(preferences: UiSettings) => void>;
  readonly #unsubscribe: ReadonlyArray<() => void>;
  #saving: Promise<void> = Promise.resolve();

  constructor(storage: PreferencesStorage, target: PageEvents, onError: (message: string) => void) {
    this.#storage = storage;
    this.#target = target;
    this.#onError = onError;
    this.#debouncer = new Debouncer(
      (preferences: UiSettings) => {
        this.#save(preferences);
      },
      { wait: SAVE_DELAY_MS },
    );
    this.#unsubscribe = [usePreferencesStore.subscribe(this.#changed), useDraftStore.subscribe(this.#changed)];
    this.#target.addEventListener("pagehide", this.#flush);
  }

  dispose(): void {
    this.#target.removeEventListener("pagehide", this.#flush);
    this.#unsubscribe.forEach((unsubscribe) => {
      unsubscribe();
    });
    this.#flush();
  }

  readonly #changed = (): void => {
    this.#debouncer.maybeExecute(currentPreferences());
  };

  readonly #flush = (): void => {
    this.#debouncer.flush();
  };

  #save(preferences: UiSettings): void {
    this.#saving = this.#saving.then(async () => this.#store(preferences));
  }

  async #store(preferences: UiSettings): Promise<void> {
    try {
      await this.#storage.save(preferences);
    } catch (error: unknown) {
      this.#onError(errorMessage(error));
    }
  }
}

export interface PageEvents {
  addEventListener(type: "pagehide", listener: () => void): void;
  removeEventListener(type: "pagehide", listener: () => void): void;
}
