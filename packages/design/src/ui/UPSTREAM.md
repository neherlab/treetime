# Vendored shadcn/ui components

Most files in this directory started as [shadcn/ui](https://ui.shadcn.com) components. The design package owns them: they follow the project lint and format rules like any other source file. The components that the apps already use keep only the parts the apps use; the others keep the registry API, because the package offers them to every site and tool.

The package wrote these components itself, without a registry source: `copy-button.tsx`, `data-table.tsx` with `sort-descriptor.ts`, `description-list.tsx`, `meter.tsx` (Base UI `Meter`, in the style of `progress.tsx`), `stat.tsx` (a value with its uncertainty and label), and `wordmark.tsx` (a name in the monospace face, optionally with `@` and a place in the style of `neherlab@biozentrum`).

## Source

- Registry: `https://ui.shadcn.com/r/styles/base-vega/<name>.json`, the `files[0].content` field, fetched on 2026-09-28 (`dropdown-menu` on 2026-10-24; `avatar`, `breadcrumb`, `hover-card`, `pagination`, `popover`, `scroll-area`, `select`, and `skeleton` on 2026-10-29)
- Style: `base-vega`, built on Base UI (`@base-ui/react`)
- Stylesheet: `shadcn.css` is `dist/tailwind.css` of the `shadcn` npm package 4.21.0, copied verbatim
- License: MIT (<https://github.com/shadcn-ui/ui/blob/main/LICENSE.md>)

## Rewrites applied to every fetched file

These are the rewrites the shadcn CLI applies when it installs a component:

- **Imports**: `@/registry/base-vega/ui/<name>` becomes `./<name>`, and `cn` comes from `./cn`
- **Icons**: each `IconPlaceholder` element becomes the Iconify Lucide icon named in its `lucide` attribute, imported as `import <Name>Icon from "~icons/lucide/<kebab-name>"` with `aria-hidden`, for example `lucide="ChevronDownIcon"` becomes `import ChevronDownIcon from "~icons/lucide/chevron-down"`
- **Classes**: `cn-font-heading` becomes `font-heading`, and the other `cn-*` marker classes are removed
- **Comments**: removed, as in all TypeScript source

## Local changes

A component refreshed from the registry needs these changes again.

- **All files**: exports and variants the app does not use are removed; `knip` reports any that become unused
- **`alert.tsx`**: a `warning` variant with a warning-colored border, background and icon, for notices that let the run go on, and a `success` variant in the success colors. The alert names its variant in `data-variant`, and the `role` prop sets the role, `alert` by default; static notes such as docs callouts pass `note`, so screen readers do not announce them as live alerts
- **`card.tsx`**: the card has the theme corner radius, a border-colored ring and no shadow
- **`accordion.tsx`**: the panel uses a 100 ms height transition instead of the keyframe animation, so a panel that mounts open appears without animation
- **`chart.tsx`**: no chart config, generated style element or legend; charts take their colors from `runs/palette.ts`. The tooltip frame is exported for custom tooltips, and the tooltip marks each series with an SVG swatch instead of inline styles
- **`command.tsx`**: built on Base UI Autocomplete instead of `cmdk`, with the class names of the registry component. `Command` is an always-open inline list that highlights the first match
- **`sidebar.tsx`**: only the left, off-canvas sidebar with the mobile sheet. No rail, icon mode, floating or inset variants, sub-menus, menu actions, badges, skeletons, tooltips or cookie persistence. The desktop sidebar has a resize handle on its right edge (drag, arrow keys, Home and End; double-click resets). Its width, 340 px by default and bounded by `sidebar-width.ts`, comes from the `width` and `onWidthChange` props, which the app keeps with its preferences, and is set as the only inline custom property; the mobile sheet width is a Tailwind class. `useMediaQuery` from `@mantine/hooks` detects the mobile layout, and `useHotkey` from `@tanstack/react-hotkeys` binds Mod+B
- **`toggle-group.tsx`**: holds the toggle variants (`toggle.tsx` is not vendored), sets the gap with a class for each supported spacing (0, 1, 2) instead of an inline custom property, and memoizes its context value
- **`input-group.tsx`**: the group is a `label` element, so a click anywhere in it focuses the input; the addon has no click handler and neither element has `role="group"`
- **`field.tsx`**: the field has no `role="group"`, which gave assistive technology an unnamed group
- **`toast.tsx`**: the default action and close buttons are module constants
- **`dropdown-menu.tsx`**: only the menu, its trigger, its content and plain items, without the `inset` and `variant` props of the items
- **`label.tsx`, `spinner.tsx`, `input-group.tsx`**: each carries one lint suppression with its reason
- **`breadcrumb.tsx`**: the current page is plain text with `aria-current="page"`, without `role="link"` and `aria-disabled`, because it is not a link
- **`breadcrumb.tsx`, `pagination.tsx`**: the ellipsis icon is Lucide `ellipsis`, the current name of `more-horizontal`; `pagination.tsx` compares the optional `isActive` with `true`
