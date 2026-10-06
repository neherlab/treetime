# Chart axis labels are passed as `label` props

Recharts axes (`XAxis`, `YAxis`) receive their title through the `label` prop, built by `bottomAxisLabel()` and `leftAxisLabel()` in `packages/app-ui/src/runs/palette.ts`, never as a `<Label>` child element.

## Reason

Recharts 3 wraps `XAxis` and `YAxis` in `memo` with a comparison that checks `label` shallowly but `children` by identity. A `<Label>` child is a new element on every render of the chart, so every chart re-render also re-renders both axes. An axis that re-renders builds new settings and replaces them in the chart store, which recomputes the scales and the bar geometry. `AnimatedItems` keys its items by an animation id that changes with the geometry array, so React remounts every bar, also with `isAnimationActive={false}`.

On the homoplasy page of `data/sc2/4500`, the genome chart re-renders on every tree recolor to move its selection line. With `<Label>` children, each recolor remounted its 4,719 bars (about 46,000 React fibers), which took about 1.5 s of a 2.5 s recolor in the dev build.

The label props objects hold only primitive values, so the shallow comparison finds two builds of one label equal, also when the label text depends on a prop.

## Implementation

- `packages/app-ui/src/runs/palette.ts`: `bottomAxisLabel()`, `leftAxisLabel()`
- Charts: `GenomeSitesChart.tsx`, `SiteHitsChart.tsx`, `RootToTipPlot.tsx`, `ShiftPlot.tsx`, `SkylinePlot.tsx` in `packages/app-ui/src/runs/`
