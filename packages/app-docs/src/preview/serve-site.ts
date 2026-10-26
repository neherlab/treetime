import { fileURLToPath } from "node:url";

import { locales } from "@neherlab/app-ui/i18n";

import { siteFile, type SiteFolders } from "./site-files";

const PORT = Number(process.env["TREETIME_DOCS_SERVE_PORT"] ?? "6700");

const FOLDERS: SiteFolders = {
  base: `/${process.env["TREETIME_DOCS_VERSION"] ?? "dev"}`,
  build: fileURLToPath(new URL("../../dist", import.meta.url)),
  root: fileURLToPath(new URL("../../dist-root", import.meta.url)),
  locales,
};

const server = Bun.serve({
  port: PORT,
  hostname: "127.0.0.1",
  async fetch(request) {
    const target = siteFile(new URL(request.url).pathname, FOLDERS);

    if (target.kind === "file") {
      for (const candidate of target.candidates) {
        const file = Bun.file(candidate);

        if (await file.exists()) {
          return new Response(file, { headers: { "Cache-Control": "no-cache" } });
        }
      }
    }

    const page = target.kind === "missing" ? target.page : target.notFound;

    return page === undefined
      ? new Response("Not found", { status: 404 })
      : new Response(Bun.file(page), { status: 404 });
  },
});

process.stdout.write(`Serving the built docs site at ${server.url.href.replace(/\/$/u, "")}${FOLDERS.base}\n`);
