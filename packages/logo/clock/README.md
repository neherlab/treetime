# Clock logo

A clock with a phylogenetic tree on its dial: the tree runs from its root at 9 o'clock to its tips at 1 to 5 o'clock, and amber needles above it point to 10 and 12 o'clock. A static version and an animated loading indicator.

## Regenerate

```bash
packages/logo/clock/generate
```

`generate` runs the whole pipeline:

1. `tools/trace`: `input/logo-concept-02.png` to the mask `output/logo-concept-02-mask.svg`. Skipped when the PNG is absent, and the committed mask is used
2. `compose`: the mask to the logos in `output/`, copied to `packages/app-ui/public/`
3. `tools/icons`: the logos to the favicons, the app icons, and the web app manifest in `packages/app-ui/public/`
4. `animate`: the loading indicator to `output/logo-loading.gif`

It runs on the host when ImageMagick (`convert`), gifsicle, potrace, `rsvg-convert`, SVGO, and `xmllint` are installed, and in the Docker image of `tools/` otherwise. `--docker` forces the image; `--host` forces the host tools.

## Files

`input/` holds the hand-made sources, `output/` everything the pipeline generates. Edit inputs and scripts, never outputs: `generate` overwrites them.

- `generate`: runs the whole pipeline
- `compose`, `animate`: the pipeline steps of this design
- `loading.html`: page with the centered loading indicator, for a quick look in a browser
- `input/logo-concept-01.png`, `input/logo-concept-02.png`: logo concepts, colored shape on a light background
- `output/logo-concept-02-mask.svg`: traced mask of `logo-concept-02.png`, a single black path in `viewBox` coordinates, without transforms
- `output/logo.svg`: the logo with full detail, for sizes from 64 px
- `output/logo-small.svg`: the logo drawn for 16 to 48 px, for the app header and the favicons
- `output/logo-loading.svg`: the animated loading indicator
- `output/logo-loading.gif`: the loading indicator as an animated GIF, for places without SVG support, such as chat messages and issue comments

## Compose

`compose <mask.svg> <output directory>` draws the clock:

- **Dial**: deep navy, with a rim in the tree gradient, so the logo keeps its outline on light and dark backgrounds
- **Tree**: gradient from cyan at the root to violet at the tips, as time runs from the past to the present. The root tip sits at the height of the 9 o'clock tick; the branches fan out toward the 1 to 5 o'clock ticks
- **Needles**: amber, at 10 and 12 o'clock, above the tree, with an outline in the dial color, like hands over a printed dial plate
- **Ticks**: 12, longer and wider at 12, 3, 6, and 9 o'clock

Details thinner than one screen pixel blur into grey, jagged edges, so `logo-small.svg` and `logo-loading.svg` are drawn for small sizes: only the 12, 3, 6, and 9 o'clock ticks, a wider and brighter rim, a tree thickened by an outline in its own gradient, and wider needles. At 24 px, the full-detail ticks would be less than half a pixel wide.

`logo-loading.svg` turns the minute needle once per 1.2 s and the hour needle 12 times slower, with CSS inside the SVG. The animation also runs when the SVG is loaded with `<img>`, and it stops when the system requests reduced motion.

`compose` measures the bounding box and the root tip of the mask, so the placement follows a re-traced mask. Ticks and needles are lines computed from their angles, and SVGO applies the tree placement to the path coordinates in a separate pass, so the outputs contain no transforms. The colors, dimensions, and detail levels are variables at the top of `compose`.

## Animate

`animate <logo-loading.svg> <output.gif>` renders one full turn of the hour needle, so the GIF loops without a jump: 14.4 s at 20 frames per second, 160 px, on a transparent background. Each frame replaces the CSS animation with the needle angles of its time as `rotate` attributes, because `rsvg-convert` draws a single frame and ignores `transform-origin`. gifsicle reduces all frames to one shared palette and stores only the pixels that change between frames. `SIZE` and `FPS` override the size and frame rate.

ImageMagick 6 shows needle trails when it reassembles the frames of this GIF (`convert -coalesce`); browsers and gifsicle play it correctly.
