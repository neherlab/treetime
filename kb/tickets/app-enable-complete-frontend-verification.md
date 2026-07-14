# Enable complete frontend verification

Make every frontend package’s declared checks exercise the code that ships.

## Required changes

- Include the desktop renderer TypeScript configuration in `typecheck`.
- Add Vite/CSS ambient declarations without weakening strict typing.
- Add meaningful web bridge, request construction, cancellation, and error tests.
- Make formatting checks pass and keep them in the workspace verification command.

## Validation

- Typecheck contracts, UI, desktop main process, desktop renderer, and web.
- Run frontend lint, format check, and Vitest with collected tests.
- Full Rust checks for shared contract generation.

## Related issues

- Source: [kb/issues/M-app-frontend-verification-excludes-renderer.md](../issues/M-app-frontend-verification-excludes-renderer.md)
