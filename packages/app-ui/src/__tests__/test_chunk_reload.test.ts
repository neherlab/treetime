import { type ChunkErrorTab, reloadOnChunkError, type ReloadTarget } from "../chunk-reload";

function page(stored: string | null) {
  const listeners: Array<(event: Event) => void> = [];
  const storage = new Map<string, string>(stored === null ? [] : [["treetime:chunk-reload-at", stored]]);
  let reloads = 0;

  const target: ReloadTarget = {
    addEventListener: (_type, listener) => {
      listeners.push(listener);
    },
    sessionStorage: {
      getItem: (key) => storage.get(key) ?? null,
      setItem: (key, value) => {
        storage.set(key, value);
      },
    },
    location: {
      reload: () => {
        reloads += 1;
      },
    },
  };

  const fail = (): boolean => {
    const event = new Event("vite:preloadError", { cancelable: true });
    listeners.forEach((listener) => listener(event));

    return event.defaultPrevented;
  };

  return { target, fail, reloads: () => reloads, stored: () => storage.get("treetime:chunk-reload-at") };
}

function tab(holdsWork: boolean) {
  let notices = 0;

  const tab: ChunkErrorTab = {
    holdsWork: () => holdsWork,
    showReloadNotice: () => {
      notices += 1;
    },
  };

  return { tab, notices: () => notices };
}

describe("chunk_reload", () => {
  test.each([
    { name: "first failure", stored: null, now: 50_000, expected: { prevented: true, reloads: 1, stored: "50000" } },
    {
      name: "failure after the interval",
      stored: "40000",
      now: 50_000,
      expected: { prevented: true, reloads: 1, stored: "50000" },
    },
    {
      name: "failure within the interval",
      stored: "45000",
      now: 50_000,
      expected: { prevented: false, reloads: 0, stored: "45000" },
    },
  ])("$name", ({ stored, now, expected }) => {
    const { target, fail, reloads, stored: storedAt } = page(stored);
    reloadOnChunkError(target, tab(false).tab, () => now);

    const prevented = fail();

    expect({ prevented, reloads: reloads(), stored: storedAt() }).toStrictEqual(expected);
  });

  test("a tab that holds work shows the reload notice and does not reload by itself", () => {
    const { target, fail, reloads, stored } = page(null);
    const { tab: busy, notices } = tab(true);
    reloadOnChunkError(target, busy, () => 50_000);

    const prevented = fail();

    expect({ prevented, reloads: reloads(), stored: stored(), notices: notices() }).toStrictEqual({
      prevented: false,
      reloads: 0,
      stored: undefined,
      notices: 1,
    });
  });
});
