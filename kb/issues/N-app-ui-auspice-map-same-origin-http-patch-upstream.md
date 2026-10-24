# Auspice map upgrades the basemap URLs of an http page to https

## Problem

The map of Auspice 3.1.0 rewrites every `http://` URL of its basemap style to `https://` before MapLibre requests it, so that tile URLs of providers such as Mapbox, whose TileJSON lists them over http, load on an https page. The TreeTime apps serve their own Natural Earth basemap from the page's origin (`mapTiles` extension in `packages/app-ui/src/auspice/extensions.ts`). On a page served over http, such as the dev server or `just serve` on localhost, the rewrite sends the requests to `https://localhost:<port>`, which no server answers, and the map shows no land, borders or labels.

TreeTime fixes this in the Bun patch of Auspice, `patches/auspice@3.1.0.patch`: a URL of the page's own origin keeps its scheme, and other http URLs are still upgraded.

## Code

- Request transform: [map/map.js#L552-L560](https://github.com/nextstrain/auspice/blob/v3.1.0/src/components/map/map.js#L552-L560)

## Fix

Contribute the change to [nextstrain/auspice](https://github.com/nextstrain/auspice), then remove it from the patch after upgrading to an Auspice release that contains it. [N-app-ui-auspice-tree-unmount-patch-upstream.md](N-app-ui-auspice-tree-unmount-patch-upstream.md) describes how to edit the patch.
