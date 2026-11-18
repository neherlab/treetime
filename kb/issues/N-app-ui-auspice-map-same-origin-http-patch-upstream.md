# Auspice map upgrades the basemap URLs of an http page to https

## Problem

The map of Auspice 3.2.0 rewrites every `http://` URL of its basemap style to `https://` before MapLibre requests it, so that tile URLs of providers such as Mapbox, whose TileJSON lists them over http, load on an https page. The TreeTime apps serve their own Natural Earth basemap from the page's origin (`mapTiles` extension in `packages/app-ui/src/auspice/extensions.ts`). On a page served over http, such as the dev server or `just serve` on localhost, the rewrite sends the requests to `https://localhost:<port>`, which no server answers, and the map shows no land, borders or labels.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.2.0.patch`: a URL of the page's own origin keeps its scheme, and other http URLs are still upgraded.

## Code

- Request transform: [map/map.js#L551-L559](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/components/map/map.js#L551-L559)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.

## Proposed implementation plan

1. Check upstream state: on the current default branch of nextstrain/auspice (clone in `.repos/auspice`), check whether `transformRequest` in `src/components/map/map.js` still upgrades every `http://` URL; the local clone is at `v3.1.0`, so fetch the default branch first
2. Take the `transformRequest` hunk of `patches/auspice@3.2.0.patch` (`isOwnOrigin = url.startsWith(`${window.location.origin}/`)`, upgrade only when not own origin) as the upstream change
3. Open a pull request at nextstrain/auspice with a reproduction: a style JSON served by the page's own http origin (for example `auspice view` on `http://localhost:4000` with a local basemap) loads no tiles before the change and loads them after it, while a Mapbox style on an https page still gets its http tile URLs upgraded
4. After an Auspice release contains the change: bump `auspice` in `packages/app-ui/package.json` and the `patchedDependencies` key in the root `package.json`, regenerate the patch without this hunk (procedure in [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md)), and delete this issue

Tests: none in TreeTime, because the change lives in Auspice. The oracle is the same-origin rule of the URL standard: a URL whose origin equals `window.location.origin` loads with the page's scheme.

Validation: `just typecheck`, `just test-ts`, and a manual check that the map in the web2 dev server over http shows land, borders, and labels.
