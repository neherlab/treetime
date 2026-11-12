# The width of the options panel is not saved

The user can drag the border between the options panel and the main area of the workspace, but the width is not saved: every new page, reload, or app start opens the panel at its default width of 26 rem. The width of the sidebar, by contrast, is a preference (`UiSettings.sidebar_width`) that web2 keeps in `localStorage` and desktop in its settings file.

## Evidence

- `ResizablePanel` with `defaultSize="26rem"` and no saved layout [packages/app-ui/src/analysis/WorkspaceLayout.tsx#L85-L91](../../packages/app-ui/src/analysis/WorkspaceLayout.tsx#L85-L91)
- `UiSettings` has `sidebar_width` and no field for the options panel (`packages/app-backend/src/app_settings/settings.rs`)

## Fix direction

Add a field for the width of the options panel to `UiSettings` beside `sidebar_width`, regenerate the contracts, keep it with the shared preferences of web2 (`sharedPreferences()` in `packages/app-ui/src/preferences/storage.ts`), and restore it as the default size of the panel.

## Validation

- A store test: a changed panel width survives a save and a load of the preferences
