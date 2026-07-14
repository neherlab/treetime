# Isolate unit-test fixtures from production data

Inventory unit tests that read repository-level `data/` at runtime. Move small biological inputs into test-local committed fixtures and load them at compile time, or reclassify tests whose contract genuinely requires production datasets as integration or golden-master tests.

## Acceptance criteria

- No unit test reads repository-level production `data/` at runtime.
- The repository-wide inventory covers every unit-test source, rather than only the currently identified examples.
- Small text inputs use test-local fixtures; compressed production datasets are not required by unit tests at runtime.
- Shared fixtures are stored once and reused by both affected flu tests.
- The biological scenarios and asserted properties remain unchanged.
- Temporarily making the production `data/` tree unavailable does not affect unit tests.

## Related issues

- Source: [kb/issues/N-test-filesystem-dependency.md](../issues/N-test-filesystem-dependency.md)
