import * as path from "node:path";

export const APP_SCHEME = "app";

export const APP_URL = `${APP_SCHEME}://treetime/`;

export const APP_SCHEME_PRIVILEGES = {
  standard: true,
  secure: true,
  supportFetchAPI: true,
  codeCache: true,
};

export type AppAsset = { kind: "file"; path: string } | { kind: "not-found" };

export function resolveAppAsset(url: string, root: string): AppAsset {
  const parsed = new URL(url);

  if (parsed.protocol !== `${APP_SCHEME}:` || parsed.host !== new URL(APP_URL).host) {
    return { kind: "not-found" };
  }

  const relative = decodePath(parsed.pathname);

  if (relative === undefined || relative.includes("\\")) {
    return { kind: "not-found" };
  }

  const file = path.resolve(root, relative);
  const inside = path.relative(root, file);

  if (inside === ".." || inside.startsWith(`..${path.sep}`) || path.isAbsolute(inside)) {
    return { kind: "not-found" };
  }

  if (path.extname(file) === "") {
    return { kind: "file", path: path.join(root, "index.html") };
  }

  return { kind: "file", path: file };
}

function decodePath(pathname: string): string | undefined {
  try {
    return decodeURIComponent(pathname).replace(/^\/+/u, "");
  } catch {
    return undefined;
  }
}
