import { useEffect, useMemo, useState } from "react";

import { fetchManifest, fetchVersions, type PageManifest, type VersionList } from "./versions";

class Cached<T> {
  readonly #load: () => Promise<T>;

  #value: Promise<T> | undefined;

  constructor(load: () => Promise<T>) {
    this.#load = load;
  }

  get(): Promise<T> {
    this.#value ??= this.#load();

    return this.#value;
  }
}

const versions = new Cached(fetchVersions);

const manifests = new Map<string, Cached<PageManifest | undefined>>();

export function useSiteVersions(): VersionList | undefined {
  return useLoaded(siteVersions);
}

export async function siteVersions(): Promise<VersionList | undefined> {
  return versions.get();
}

export function useManifest(version: string | undefined): PageManifest | undefined {
  const load = useMemo(() => (version === undefined ? undefined : async () => manifestOf(version)), [version]);

  return useLoaded(load);
}

export async function manifestOf(version: string): Promise<PageManifest | undefined> {
  const cached = manifests.get(version) ?? new Cached(async () => fetchManifest(version));
  manifests.set(version, cached);

  return cached.get();
}

function useLoaded<T>(load: (() => Promise<T | undefined>) | undefined): T | undefined {
  const [value, setValue] = useState<T | undefined>(undefined);

  useEffect(() => {
    let active = true;

    async function update(): Promise<void> {
      const loaded = await load?.();

      if (active) {
        setValue(loaded);
      }
    }

    void update();

    return () => {
      active = false;
    };
  }, [load]);

  return value;
}
