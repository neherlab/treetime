# Mugration tables never show state colors

`StateLabel` in `packages/app-ui/src/runs/MugrationResultsView.tsx` draws a colored dot before each state in the "State changes" and "Uncertain ancestors" tables. It takes the color from the `scale` of the Auspice coloring whose key is the mugration attribute. The Auspice file that `treetime mugration` writes has no `scale` for that coloring, so the color map is empty and no dot is drawn.

## Evidence

Mugration run on `data/dengue/1000` with `--attribute country` in the dev app: `meta.colorings` of `mugration.auspice.json` holds `{"key": "country", "type": "categorical"}` without `scale`, and neither table shows a dot. Auspice assigns the tree colors in the browser, so the tree is colored while the tables are not.

## Fix direction

Give the tables the colors that the tree uses: either write a `scale` for the attribute into the Auspice file, or read the colors Auspice assigned.

> [!IMPORTANT]
> **Decision required.** A `scale` in the Auspice file fixes the colors for every Auspice viewer of the file, while reading Auspice's colors keeps the file unchanged and couples the tables to the Auspice store. Which owner provides state colors is open.
