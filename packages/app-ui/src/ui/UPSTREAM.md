# Vendored shadcn/ui components

The files in this directory started as [shadcn/ui](https://ui.shadcn.com) components. TreeTime owns them: they follow the project lint and format rules like any other source file, and they keep only the parts the app uses.

## Source

- Registry: `https://ui.shadcn.com/r/styles/base-vega/<name>.json`, the `files[0].content` field, fetched on 2026-09-28
- Style: `base-vega`, built on Base UI (`@base-ui/react`)
- Stylesheet: `shadcn.css` is `dist/tailwind.css` of the `shadcn` npm package 4.21.0, copied verbatim
- License: MIT (<https://github.com/shadcn-ui/ui/blob/main/LICENSE.md>)

## Rewrites applied to every fetched file

These are the rewrites the shadcn CLI applies when it installs a component:

- **Imports**: `@/registry/base-vega/ui/<name>` becomes `./<name>`, and `cn` comes from `./cn`
- **Icons**: each `IconPlaceholder` element becomes the `lucide-react` icon named in its `lucide` attribute
- **Classes**: `cn-font-heading` becomes `font-heading`, and the other `cn-*` marker classes are removed
- **Comments**: removed, as in all TypeScript source

## Local changes

A component refreshed from the registry needs these changes again.

- **All files**: exports and variants the app does not use are removed; `knip` reports any that become unused
- **`card.tsx`**: the card has the theme corner radius, a border-colored ring and no shadow, and its title is bold
- **`dialog.tsx`, `sheet.tsx`, `empty.tsx`**: titles are bold, because Lato has no medium weight
- **`chart.tsx`**: no chart config, generated style element or legend; charts take their colors from `runs/palette.ts`. The tooltip frame is exported for custom tooltips, and the tooltip marks each series with an SVG swatch instead of inline styles
- **`command.tsx`**: built on Base UI Autocomplete instead of `cmdk`, with the class names of the registry component. `Command` is an always-open inline list that highlights the first match
- **`sidebar.tsx`**: only the left, off-canvas sidebar with the mobile sheet. No rail, icon mode, floating or inset variants, sub-menus, menu actions, badges, skeletons, tooltips or cookie persistence. Widths are Tailwind classes instead of inline custom properties; `useMediaQuery` from `@mantine/hooks` detects the mobile layout, and `useHotkey` from `@tanstack/react-hotkeys` binds Mod+B
- **`toggle-group.tsx`**: holds the toggle variants (`toggle.tsx` is not vendored), sets the gap with a class for each supported spacing (0, 1, 2) instead of an inline custom property, and memoizes its context value
- **`input-group.tsx`**: the group is a `label` element, so a click anywhere in it focuses the input; the addon has no click handler and neither element has `role="group"`
- **`field.tsx`**: the field has no `role="group"`, which gave assistive technology an unnamed group
- **`toast.tsx`**: the default action and close buttons are module constants
- **`label.tsx`, `spinner.tsx`, `input-group.tsx`**: each carries one lint suppression with its reason
