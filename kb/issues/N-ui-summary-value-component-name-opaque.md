# Results summary component uses the opaque `KV` name

The local React component `KV` renders one labeled summary value [`packages/app-ui/src/components/ResultsPanel.tsx#L484`](../../packages/app-ui/src/components/ResultsPanel.tsx#L484). Its call sites provide clock and pruning summary rows, but the abbreviation only hints at an implementation-level key-value pair.

## Current role

The component renders a label and one formatted result value. It is used for clock and pruning summaries beside components named for presentation roles such as `StatCard` and `ResultsPanel`. It is not a generic key-value data structure.

## Required change

Rename the component for its UI role, such as `SummaryValue`, so searches and JSX call sites reveal what it presents. Keep its semantic markup and theme-reactive styling unchanged.

## Validation

- Component tests cover the rendered label and value.
- TypeScript checking confirms every JSX reference moved with the symbol.
- UI snapshots, if present, remain unchanged apart from component names in diagnostics.
