import { join } from "node:path";

export function siteFile(urlPath: string, folders: SiteFolders): SiteFile {
  const path = decodeURIComponent(urlPath.split(/[?#]/u)[0] ?? "/");

  if (path.split("/").includes("..")) {
    return { kind: "missing", page: join(folders.build, "404.html") };
  }

  if (path === "/" || path === "/index.html") {
    return { kind: "file", candidates: [join(folders.root, "index.html")] };
  }

  if (!path.startsWith(`${folders.base}/`) && path !== folders.base) {
    return { kind: "file", candidates: [join(folders.root, path)] };
  }

  const relative = path.slice(folders.base.length).replace(/^\/+/u, "");

  return {
    kind: "file",
    candidates:
      relative === ""
        ? [join(folders.build, "index.html")]
        : [join(folders.build, relative), join(folders.build, `${relative}.html`)],
    notFound: notFoundPage(relative, folders),
  };
}

export interface SiteFolders {
  base: string;
  build: string;
  root: string;
  locales: readonly string[];
}

export type SiteFile = { kind: "file"; candidates: string[]; notFound?: string } | { kind: "missing"; page: string };

function notFoundPage(relative: string, folders: SiteFolders): string {
  const locale = folders.locales.find((candidate) => relative === candidate || relative.startsWith(`${candidate}/`));

  return join(folders.build, locale === undefined ? "404.html" : join(locale, "404.html"));
}
