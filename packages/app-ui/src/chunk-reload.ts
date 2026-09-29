import { DateTime } from "luxon";

const RELOAD_KEY = "treetime:chunk-reload-at";

const RELOAD_INTERVAL_MS = 10_000;

export interface ReloadTarget {
  addEventListener: (type: "vite:preloadError", listener: (event: Event) => void) => void;
  sessionStorage: Pick<Storage, "getItem" | "setItem">;
  location: Pick<Location, "reload">;
}

export function reloadOnChunkError(target: ReloadTarget, now: () => number = currentMillis): void {
  target.addEventListener("vite:preloadError", (event) => {
    const at = now();

    if (!canReload(target.sessionStorage.getItem(RELOAD_KEY), at)) {
      return;
    }

    event.preventDefault();
    target.sessionStorage.setItem(RELOAD_KEY, String(at));
    target.location.reload();
  });
}

function canReload(lastReload: string | null, at: number): boolean {
  return lastReload === null || at - Number(lastReload) >= RELOAD_INTERVAL_MS;
}

function currentMillis(): number {
  return DateTime.now().toMillis();
}
