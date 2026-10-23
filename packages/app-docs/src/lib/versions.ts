import * as z from "zod";

const FETCH_TIMEOUT_MS = 10_000;

const versionEntrySchema = z.object({
  version: z.string().regex(/^[a-z0-9][a-z0-9.-]*$/u),
  label: z.string().min(1),
  preferred: z.boolean().optional(),
});

export const versionListSchema = z
  .array(versionEntrySchema)
  .min(1)
  .refine((entries) => entries.filter((entry) => entry.preferred === true).length === 1, {
    message: "exactly one version is preferred",
  })
  .refine((entries) => new Set(entries.map((entry) => entry.version)).size === entries.length, {
    message: "each version appears once",
  });

export type VersionList = z.infer<typeof versionListSchema>;

export type VersionEntry = z.infer<typeof versionEntrySchema>;

export const pageManifestSchema = z.object({
  version: z.string(),
  pages: z.record(z.string(), z.array(z.string())),
});

export type PageManifest = z.infer<typeof pageManifestSchema>;

export function preferredVersion(versions: VersionList): VersionEntry {
  const preferred = versions.find((entry) => entry.preferred === true);

  if (preferred === undefined) {
    throw new Error("the version list has no preferred version");
  }

  return preferred;
}

export async function fetchVersions(): Promise<VersionList | undefined> {
  return fetchValidated("/versions.json", versionListSchema);
}

export async function fetchManifest(version: string): Promise<PageManifest | undefined> {
  return fetchValidated(`/${version}/pages.json`, pageManifestSchema);
}

async function fetchValidated<T>(url: string, schema: z.ZodType<T>): Promise<T | undefined> {
  try {
    const response = await fetch(url, { signal: AbortSignal.timeout(FETCH_TIMEOUT_MS) });

    if (!response.ok) {
      return undefined;
    }

    const parsed = schema.safeParse(await response.json());

    return parsed.success ? parsed.data : undefined;
  } catch {
    return undefined;
  }
}
