# Logo

TreeTime logo designs and the tools that build them. Each design lives in its own directory: its scripts at the top, hand-made sources in `input/`, and generated files in `output/`. The tools are shared.

- `tools/`: shared image tools and their Docker image: tracing a raster into an SVG mask, rasterizing favicons and app icons. See `tools/README.md`
- `clock/`: a clock with a phylogenetic tree on its dial, static and animated. See `clock/README.md`

The design whose `generate` ran last owns the favicons, the app icons, and the logos in `packages/app-ui/public/`, the static directory that the web app and the desktop app both serve.
