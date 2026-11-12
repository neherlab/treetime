# The log view and the setting help each define the same Markdown plugin

The Log view of a run and the shared setting help render Markdown from Rust with the same plugin list, `remarkGfm` and a local `remarkHtmlAsText()` that shows raw HTML as text. Each file defines `remarkHtmlAsText()` and the plugin list itself, so a change of the rendering rules (for example of the sanitizing or of raw HTML) reaches one renderer and not the other.

## Evidence

- `REMARK_PLUGINS` and `remarkHtmlAsText()` in [packages/app-ui/src/runs/LogMessageView.tsx#L6](../../packages/app-ui/src/runs/LogMessageView.tsx#L6) and [#L50](../../packages/app-ui/src/runs/LogMessageView.tsx#L50)
- The same in [packages/app-ui/src/settings/help.tsx#L16](../../packages/app-ui/src/settings/help.tsx#L16) and [#L165](../../packages/app-ui/src/settings/help.tsx#L165)

## Fix direction

Move the plugin lists and `remarkHtmlAsText()` into one shared Markdown module of `app-ui` (or one Markdown component), and render both views through it; add it to the utility index `docs/en/development/utilities.md`.

## Validation

- One definition of the Markdown rules; a test of the shared module shows raw HTML as text
