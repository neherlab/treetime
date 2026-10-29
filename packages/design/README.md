# @neherlab/design

The design language of the neherlab sites and tools: the look of [neherlab.org](https://neherlab.org) as tokens, themes, fonts, and React components, in a light and a dark theme. The developer guide page [`docs/en/development/design-language.mdx`](../../docs/en/development/design-language.mdx) describes the registers and the principles and shows a specimen of every token and component.

| Path                                    | Contents                                                                                       |
| --------------------------------------- | ---------------------------------------------------------------------------------------------- |
| `src/tokens/tokens.ts`                  | The token source: every color as a light and dark OKLCH pair, typefaces, corner radii, motion |
| `src/tokens/color.ts`                   | Conversion of tokens to CSS, hex, and bytes, and WCAG contrast                                |
| `src/css/tokens.css`                    | CSS custom properties of every token (generated)                                              |
| `src/css/tailwind.css`                  | The Tailwind CSS theme, variants, and base styles; imports the generated `tailwind-theme.css` |
| `src/css/docusaurus.css`                | Preset for Docusaurus (Infima variables)                                                       |
| `src/css/reveal.css`                    | Theme for reveal.js 4 and 5, loaded after a stock theme                                        |
| `src/fonts.ts`                          | Self-hosted Source Serif 4, IBM Plex Sans, and IBM Plex Mono                                   |
| `src/ui/`                               | React components on Base UI and React Aria Components; `UPSTREAM.md` lists their sources       |
| `src/specimen/`                         | The specimen of the docs page                                                                  |
| `matplotlib/neherlab.mplstyle`          | Style for matplotlib figures (generated)                                                       |
| `build/generate.ts`                     | Writes the generated files, and the chart palette of `treetime-charts`                         |

After a change of `src/tokens/tokens.ts`, run `just gen design`. `just test-ts` checks the contrast of every text, control, and chart color in both themes, and that every stylesheet of the package reads only defined tokens.
