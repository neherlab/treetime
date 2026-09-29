# Logo tools

Image tools shared by the logo designs in `packages/logo/`.

## Files

- `run`: runs a command in the Docker image of the image tools
- `Dockerfile`, `Dockerfile.dockerignore`: the image, with ImageMagick, gifsicle, potrace, SVGO, `xmllint`, and `rsvg-convert`
- `package.json`, `package-lock.json`: the pinned SVGO version and its locked dependencies
- `svgo.config.mjs`: SVGO settings: 2 decimals, and removal of the fixed width and height
- `trace`: traces a logo raster into a single-path SVG mask
- `icons`: rasterizes the favicons, the app icons, and the web app manifest from a logo

## Run in Docker

```bash
packages/logo/tools/run <command> [args...]
```

`run` builds the image on first use and again when the Dockerfile, the dockerignore file, or the npm files change. It mounts the repository at its host path, with the git metadata read-only, and keeps the current directory, so paths work the same in the container and on the host. The container runs as the current user, without network access and without Linux capabilities. `NETWORK=bridge` enables network access for one command.

## Trace

`trace <input.png> <output.svg>` traces the colored shape of a logo raster on a light background:

1. **Shape mask**: for each pixel, the maximum channel distance from the background color (the top-left pixel). Light and dark shape colors both separate from the background this way
2. **Upscale**: the raster is upscaled before thresholding, so the anti-aliased edges become sub-pixel geometry for the tracer
3. **Edge smoothing**: the binary mask is blurred and thresholded again at 50%. This removes edge noise of the raster without moving straight edges
4. **Vectorization**: potrace fits Bezier curves into one flat path
5. **Cleanup**: SVGO applies the potrace group transform to the path coordinates, removes the doctype, metadata, and fixed size, and rounds to 2 decimals

potrace computes in a y-up coordinate system with the origin at the lower-left corner and rounds points to integers of 1/`unit` pixel. Its SVG output therefore always wraps the path in `translate(0,H) scale(s,-s)`, and no potrace option avoids it. SVGO multiplies this transform into the path data, so the numbers in `d` are coordinates of the `viewBox`. With the 4x upscale, coordinates are multiples of 0.25 source pixels, so 2 decimals round nothing away.

Settings, as environment variables of `trace`:

- `SCALE`: upscale factor before thresholding. Default `4`
- `THRESHOLD`: background distance that counts as shape. Default `35%`
- `BLUR`: edge smoothing sigma in upscaled pixels. Default `6`. Lower values keep small kinks of the raster edge; higher values round sharp corners and shorten thin tips

To preview a mask as PNG:

```bash
packages/logo/tools/run rsvg-convert -w 1076 -b white packages/logo/clock/output/logo-concept-02-mask.svg -o preview.png
```

## Icons

`icons <logo.svg> <logo-small.svg> <output directory>` rasterizes the favicons and app icons: the favicons from the small logo, the larger icons from the full-detail logo. The opaque icons use the first gradient color of the logo's first circle as background, which also sets the manifest theme and background colors.

- `favicon.svg`, `favicon.ico`: the browser tab icons; the ICO holds 16, 32, and 48 px images for browsers and tools without SVG favicons
- `apple-touch-icon.png`: 180 px home screen icon for iOS, opaque
- `icon-192.png`, `icon-512.png`: web app manifest icons with a transparent background
- `icon-maskable-512.png`: web app manifest icon for Android launchers, which crop icons to a circle, a rounded square, or another shape. The logo fills the centered safe circle of 80% diameter
- `manifest.webmanifest`: web app name, colors, and icons

## Update SVGO

Change the version in `package.json`, then regenerate the lockfile with network access. `--before` applies the seven-day supply-chain quarantine to every locked package:

```bash
NETWORK=bridge packages/logo/tools/run npm install --prefix=packages/logo/tools --package-lock-only --ignore-scripts --no-audit --no-fund --cache=/tmp/npm --before="$(date -u -d '7 days ago' +%Y-%m-%dT%H:%M:%SZ)"
```

The Debian 12 image ships Node.js 18; every locked package must support it (`engines` in `package-lock.json`).
