---
title: Web app
description: Run, serve, and deploy the web app with its API server.
---

`just up` starts the API server and the Vite dev server in the foreground until Ctrl-C. `just health` and `just status` report whether they run and which commit they were started from.

`just serve` builds the web app and the API server for production and runs the server in the foreground until Ctrl-C. The server serves the built web app (`packages/app-web/dist`, through the `STATIC_DIR` environment variable) and the API from one port.

Each checkout takes its ports from its port slot, a group of 10 ports in the block 14000-14999, so worktrees run side by side. The main checkout holds slot 0 (ports 14000-14007). A linked worktree claims the lowest slot that no other checkout of the repository holds and whose ports have no listener, and keeps it until the worktree is removed; the claims are files in the git common directory (`dev/lib/ports.sh`). `just status` prints the slot and its ports. `TREETIME_API_PORT`, `TREETIME_WEB_PORT`, `TREETIME_SERVE_PORT`, and the other `TREETIME_*_PORT` variables override single ports.

A route launcher on the host can publish the web dev server through a local HTTPS proxy, such as [portless](https://github.com/vercel-labs/portless), under a hostname such as `https://treetime.localhost`. The launcher starts `just up` as its child process, assigns the hostnames, and removes the route when the command ends. The checkout and the launcher exchange two files in `tmp/portless/`:

- `web.hostnames.json`: `{"hostnames": [...], "urls": [...]}`, written by the launcher before it starts the command. `just status` prints the URLs. Without this file, the server has no public hostname
- `web.json`: `{"port": <n>}`, written by `just up` after the web server answers, through a temporary file and a rename, so the launcher never reads a partial file. `just up` writes it also when no launcher runs

The checkout never reads or writes the state of the proxy itself. The Vite server accepts every `*.localhost` hostname, so each assigned name works without configuration.

## Storage

The web app has no database. The server keeps each run in its own folder `<runs-dir>/<run-id>/`: the run record `run.json`, the event log `events.jsonl`, the uploaded files in `inputs/`, and the results in `out/`. The browser keeps the sidebar width in `localStorage` under the key `treetime-preferences`, so every tab shares it, and the unfinished form in `sessionStorage` under the key `treetime-draft`, so each tab keeps its own form through a reload. The theme stays with `next-themes`, which stores it in `localStorage` under the key `theme` and applies a change in one tab to the other tabs at once. The server reads run inputs only from the `inputs/` folders of the runs; the web app uploads a chosen example file into the run like a user file (see [Example data](/development/desktop-app#example-data)).

Each app mode keeps its runs in its own directory:

- `tmp/app/web-dev/runs`: development servers (`just up`)
- `tmp/app/web-prod/runs`: production server (`just serve`)
- `tmp/app/treetime-dev/runs`: desktop app in development mode (`just desktop`); see [Desktop app](/development/desktop-app) for its other folders
- `tmp/app/treetime-prod/runs`: production build of the desktop app started from the checkout (`just desktop-prod`)
- `<app folder>/runs`: installed desktop app, unless the settings or `TREETIME_RUNS_DIR` name another folder; see [Desktop app](/development/desktop-app)

A deployment passes `--runs-dir` to `treetime-server` itself.

## Event streams and HTTP/2

The web app follows changes through server-sent event streams: `GET /api/events` for changes to any run, and `GET /api/runs/{id}/events` for each run it follows. Every open stream holds one connection. Over HTTP/1.1, browsers allow at most 6 connections per host, shared by all tabs, so more open streams stall every further request to the server. Over HTTP/2, all requests share one connection.

The API server speaks HTTP/1.1 without TLS, and browsers use HTTP/2 only over TLS. A deployment therefore puts a reverse proxy in front of the server that terminates TLS and serves HTTP/2 to the browser, for example nginx with `listen 443 ssl; http2 on;`. In development, the HTTPS proxy of a route launcher (see above) does this. Without such a proxy, at most 6 streams stay open per host across all tabs.

## Caching, compression, and allowed hosts

`treetime-server` serves the UI build (`STATIC_DIR`) and the API with these rules, so a proxy in front of it needs no caching or compression settings of its own:

- **UI files**: the file names in `assets/` contain a hash of their content, so they are cached for a year (`Cache-Control: public, max-age=31536000, immutable`). Every other file, `index.html` included, is revalidated on each use (`no-cache`), so a new deploy reaches the browser at the next page load. Only page navigations (`Accept: text/html`) fall back to `index.html`; a missing asset or API path answers 404
- **Compression of UI files**: `bun run build:web` writes a gzip (level 9) and a brotli (quality 9) copy next to each text file of the build, and the server sends the copy the browser accepts, with `Vary: Accept-Encoding`. Brotli quality 11 would make the files about 8 % smaller but adds seconds to every build
- **Compression of API responses**: zstd, or gzip for browsers without zstd, at the fastest level. Run outputs of large alignments shrink far more with zstd than with gzip, because gzip only finds repeats within 32 KB, less than one genome. Event streams and already compressed files (`.xz`, `.gz`, `.zst`, fonts) stay uncompressed
- **API caching**: API responses are not stored (`no-store`). The results, the Auspice document, and the comparison of finished runs carry an `ETag`, so a browser that asks again gets `304 Not Modified` without the server rebuilding the answer. The tag changes when the server restarts
- **Downloads**: output files and run archives stream to the browser, so the server holds only a small buffer per download. Archives store their files uncompressed and rely on the HTTP compression above
- **Allowed hosts**: the server answers only requests for `localhost`, `*.localhost`, `127.0.0.1`, `[::1]`, and the hosts given with `--allowed-host` (repeatable) or `ALLOWED_HOSTS` (comma-separated); other hosts get HTTP 403. This stops other websites from reaching a local server through DNS rebinding. The server sends no CORS headers, because the UI calls it from its own origin
