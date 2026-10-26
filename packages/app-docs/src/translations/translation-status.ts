import { createHash } from "node:crypto";

import { baseLocale, locales } from "@neherlab/app-ui/i18n";
import * as z from "zod";

const TRANSLATED_FILE = /\.(?:md|mdx|json)$/u;

const FRONT_MATTER = /^---\r?\n(?<yaml>[\s\S]*?)\r?\n---(?:\r?\n|$)/u;

const sourceHashSchema = z.looseObject({ sourceHash: z.string().optional() });

export function translationStatus(files: ReadonlyMap<string, string>): TranslationEntry[] {
  const sources = [...files.keys()].filter((path) => path.startsWith(`${baseLocale}/`) && TRANSLATED_FILE.test(path));
  const targetLocales = locales.filter((locale) => locale !== baseLocale);

  const pending = sources.flatMap((source) => {
    const page = source.slice(baseLocale.length + 1);
    const hash = contentHash(files.get(source) ?? "");

    return targetLocales.flatMap((locale): TranslationEntry[] => {
      const target = `${locale}/${page}`;
      const content = files.get(target);

      if (content === undefined) {
        return [{ status: "missing", path: target, sourceHash: hash }];
      }

      return recordedHash(target, content) === hash ? [] : [{ status: "stale", path: target, sourceHash: hash }];
    });
  });

  const orphans = [...files.keys()].flatMap((path): TranslationEntry[] =>
    isOrphan(path, files, targetLocales) ? [{ status: "orphan", path, sourceHash: undefined }] : [],
  );

  return [...pending, ...orphans].toSorted((left, right) => left.path.localeCompare(right.path));
}

export function contentHash(content: string): string {
  return createHash("sha256").update(content).digest("hex").slice(0, 16);
}

export function formatEntry({ status, path, sourceHash }: TranslationEntry): string {
  return sourceHash === undefined ? `${status} ${path}` : `${status} ${path} ${sourceHash}`;
}

export interface TranslationEntry {
  status: "missing" | "stale" | "orphan";
  path: string;
  sourceHash: string | undefined;
}

function isOrphan(path: string, files: ReadonlyMap<string, string>, targetLocales: readonly string[]): boolean {
  const translated = TRANSLATED_FILE.test(path) && targetLocales.some((locale) => path.startsWith(`${locale}/`));

  return translated && !files.has(`${baseLocale}/${path.slice(path.indexOf("/") + 1)}`);
}

function recordedHash(path: string, content: string): string | undefined {
  const data: unknown = path.endsWith(".json")
    ? JSON.parse(content)
    : Bun.YAML.parse(FRONT_MATTER.exec(content)?.groups?.["yaml"] ?? "");

  return sourceHashSchema.safeParse(data).data?.sourceHash;
}
