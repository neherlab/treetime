# App supports nested settings but no command has one

## Problem

The app's setting machinery handles nested settings: key paths such as `group.field` in the catalog, nested YAML in
the config code view, nested comparisons in the run settings and flag lookup of nested schema properties
([packages/app-backend/src/config/](../../packages/app-backend/src/config/),
[packages/app-ui/src/settings/](../../packages/app-ui/src/settings/)). Its only users were the `clock` groups
`branch_split` and `clock_regression`, removed with the retired clock flags
([clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md)). No command
schema has a nested setting now, and the tests that exercised the nested paths with those groups were dropped or
moved to flat settings (commit "fix(app): drop the covariation setting and the nested-setting examples it took with
it").

> [!IMPORTANT]
> **Decision required.** Keep the nested-setting support, tested on a synthetic schema or config, or remove it
> until a command needs it again. Keeping it leaves untested code paths in the backend, contracts and UI; removing
> it simplifies them and must be redone if a nested setting returns.

## Validation

Either the nested paths are tested again without a real nested setting, or the code that handles them is gone and
the remaining tests pass.
