import { mkdir, readdir, readFile, rm, writeFile } from "node:fs/promises";
import { join, relative } from "node:path";
import { fileURLToPath } from "node:url";

import { baseLocale } from "@neherlab/app-ui/i18n";
import type { AstroIntegration } from "astro";

import { redirectPage } from "../lib/redirect-page";
import { docsPath } from "../lib/urls";
import type { PageManifest, VersionList } from "../lib/versions";
import { preferredVersion } from "../lib/versions";
import { isContentPage, manifestPages } from "./built-pages";
import { writeCompressedCopies } from "./compress";
import { cspViolations } from "./csp-check";
import { writeSearchIndex } from "./search-index";

const SEARCH_FOLDER = "pagefind";

export function postBuild(options: PostBuildOptions): AstroIntegration {
  return {
    name: "treetime-docs-post-build",
    hooks: {
      "astro:build:done": async ({ dir, logger }) => {
        const root = fileURLToPath(dir);
        const pages = await builtPages(root);

        await writeSearchIndex(root, options.base, pages.filter(isContentPage));
        await writeManifest(root, options.version, pages);
        await writeSiteRoot(fileURLToPath(new URL("../dist-root/", dir)), options);
        const compressed = await writeCompressedCopies(root, [SEARCH_FOLDER]);
        await checkContentSecurityPolicy(root, pages);

        logger.info(`indexed ${pages.filter(isContentPage).length} pages and compressed ${compressed} files`);
      },
    },
  };
}

export interface PostBuildOptions {
  base: string;
  version: string;
  versions: VersionList;
  site: string;
}

async function builtPages(root: string): Promise<string[]> {
  const entries = await readdir(root, { recursive: true, withFileTypes: true });

  return entries
    .filter((entry) => entry.isFile() && entry.name.endsWith(".html"))
    .map((entry) => relative(root, join(entry.parentPath, entry.name)))
    .filter((file) => !file.startsWith(`${SEARCH_FOLDER}/`))
    .toSorted();
}

async function writeManifest(root: string, version: string, pages: readonly string[]): Promise<void> {
  const manifest: PageManifest = { version, pages: manifestPages(pages) };

  await writeFile(join(root, "pages.json"), `${JSON.stringify(manifest)}\n`);
}

async function writeSiteRoot(folder: string, options: PostBuildOptions): Promise<void> {
  await rm(folder, { recursive: true, force: true });
  await mkdir(folder, { recursive: true });
  await writeFile(join(folder, "versions.json"), `${JSON.stringify(options.versions)}\n`);
  await writeFile(
    join(folder, "index.html"),
    redirectPage(options.site, docsPath(`/${preferredVersion(options.versions).version}`, [], baseLocale)),
  );
}

async function checkContentSecurityPolicy(root: string, pages: readonly string[]): Promise<void> {
  const violations = (
    await Promise.all(pages.map(async (page) => cspViolations(page, await readFile(join(root, page), "utf8"))))
  ).flat();

  if (violations.length > 0) {
    throw new Error(`The content security policy check failed:\n${violations.join("\n")}`);
  }
}
