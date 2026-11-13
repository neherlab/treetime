# The desktop app does not start when its settings file does not parse

The desktop back end reads `settings.yaml` (or `settings.json`) when it opens, and a file that does not parse stops the start. The settings types deny unknown fields, and the file also holds the unfinished form (`ui.draft`), whose schema changes with the form. So a change of the draft schema, or of any settings type, stops the desktop app on a machine whose settings file an older version wrote, and the user sees no app. web2 drops an unreadable draft from `sessionStorage` and starts.

## Evidence

- `AppSettingsStore::read_file()` returns the parse error with the file path [packages/app-backend/src/app_settings/store.rs#L64-L72](../../packages/app-backend/src/app_settings/store.rs#L64-L72)
- `DesktopService::open()` propagates it: `let settings = store.read()?;` [packages/app-napi/src/backend.rs#L41](../../packages/app-napi/src/backend.rs#L41)
- `AppSettings` and its parts derive with `deny_unknown_fields` [packages/app-backend/src/app_settings/settings.rs#L15-L16](../../packages/app-backend/src/app_settings/settings.rs#L15-L16)
- web2: `BrowserPreferencesStorage` parses the stored draft and reports and drops a value that does not parse (`packages/app-ui/src/preferences/storage.ts`)

## Fix direction

Read the draft apart from the folders and the preferences: a draft that does not parse is dropped with a warning in the log, as web2 does, while an error in `paths` keeps stopping the start with the file path, because the folders decide where the runs are. Test: a settings file with a draft of an unknown shape opens the app with an empty form.

## Validation

- The desktop back end opens with a settings file whose draft does not parse, and keeps its folders and theme
