# Logo

TreeTime logo, its concepts, and the tools that build it: `trace` turns a logo raster into a monochrome SVG mask, `compose` draws the clock logo and its loading indicator from the mask.

## Files

- `logo-concept-01.png`, `logo-concept-02.png`: logo concepts, colored shape on a light background
- `logo.svg`, `logo-loading.svg`: the logo, a clock with the traced tree on its dial, static and animated
- `loading.html`: page with the centered loading indicator
- `logo-concept-02-mask.svg`: traced mask of `logo-concept-02.png`, a single black path in `viewBox` coordinates, without transforms
- `trace`: tracing script
- `compose`: builds `logo.svg` and `logo-loading.svg` from a mask
- `svgo.config.mjs`: SVGO settings of `trace` and `compose`
- `run`: runs a command in the Docker image of the image tools
- `Dockerfile`, `Dockerfile.dockerignore`: the image, with ImageMagick, potrace, SVGO, `xmllint`, and `rsvg-convert`
- `package.json`, `package-lock.json`: the pinned SVGO version and its locked dependencies

## Trace a logo

```bash
packages/logo/run ./trace logo-concept-02.png logo-concept-02-mask.svg
```

`run` builds the image on first use and again when the Dockerfile, the dockerignore file, or the npm files change. It mounts this directory at `/workdir`, so paths are relative to this directory. The container runs as the current user, without network access and without Linux capabilities. `NETWORK=bridge` enables network access for one command.

To preview a mask as PNG:

```bash
packages/logo/run rsvg-convert -w 1076 -b white logo-concept-02-mask.svg -o preview.png
```

## How tracing works

1. **Shape mask**: for each pixel, the maximum channel distance from the background color (the top-left pixel). Light and dark shape colors both separate from the background this way
2. **Upscale**: the raster is upscaled before thresholding, so the anti-aliased edges become sub-pixel geometry for the tracer
3. **Edge smoothing**: the binary mask is blurred and thresholded again at 50%. This removes edge noise of the raster without moving straight edges
4. **Vectorization**: potrace fits Bezier curves into one flat path
5. **Cleanup**: SVGO applies the potrace group transform to the path coordinates, removes the doctype, metadata, and fixed size, and rounds to 2 decimals

potrace computes in a y-up coordinate system with the origin at the lower-left corner and rounds points to integers of 1/`unit` pixel. Its SVG output therefore always wraps the path in `translate(0,H) scale(s,-s)`, and no potrace option avoids it. SVGO multiplies this transform into the path data, so the numbers in `d` are coordinates of the `viewBox`. With the 4x upscale, coordinates are multiples of 0.25 source pixels, so 2 decimals round nothing away

## Settings

Set as environment variables of `trace`, for example `packages/logo/run env BLUR=4 ./trace in.png out.svg`:

- `SCALE`: upscale factor before thresholding. Default `4`
- `THRESHOLD`: background distance that counts as shape. Default `35%`
- `BLUR`: edge smoothing sigma in upscaled pixels. Default `6`. Lower values keep small kinks of the raster edge; higher values round the sharp fork corners and shorten the branch tips

## Compose the logo

```bash
packages/logo/run ./compose logo-concept-02-mask.svg logo.svg logo-loading.svg
```

`compose` draws a clock with the traced tree on its dial and writes two files:

- `logo.svg`: the static logo, needles at 10 and 12 o'clock
- `logo-loading.svg`: the loading indicator. CSS inside the SVG turns the minute needle once per 1.2 s and the hour needle 12 times slower. The animation also runs when the SVG is loaded with `<img>`, and it stops when the system requests reduced motion

Design:

- **Dial**: deep navy, with a thin rim in the tree gradient, so the logo keeps its outline on light and dark backgrounds
- **Tree**: gradient from cyan at the root to violet at the tips, as time runs from the past to the present. The root tip sits at the height of the 9 o'clock tick; the branches fan out toward the 1 to 5 o'clock ticks
- **Needles**: amber, above the tree, with an outline in the dial color, like hands over a printed dial plate
- **Ticks**: 12, longer and wider at 12, 3, 6, and 9 o'clock

`compose` measures the bounding box and the root tip of the mask, so the placement follows a re-traced mask. Ticks and needles are lines computed from their angles, and SVGO applies the tree placement to the path coordinates in a separate pass, so the output contains no transforms. The colors and dimensions are variables at the top of `compose`.

`loading.html` shows the loading indicator centered on a page that follows the light or dark color scheme of the system.

## Update SVGO

Change the version in `package.json`, then regenerate the lockfile with network access. `--before` applies the seven-day supply-chain quarantine to every locked package:

```bash
NETWORK=bridge packages/logo/run bash -c 'npm install --package-lock-only --ignore-scripts --no-audit --no-fund --cache=/tmp/npm --before="$(date -u -d "7 days ago" +%Y-%m-%dT%H:%M:%SZ)"'
```

The Debian 12 image ships Node.js 18; every locked package must support it (`engines` in `package-lock.json`).
