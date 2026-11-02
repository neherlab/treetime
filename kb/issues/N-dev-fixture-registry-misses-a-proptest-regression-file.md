# The fixture registry misses a proptest regression file

`just fixtures-check` fails on `rust`:

```text
tracked fixtures that no registry entry claims (register them in dev/registry/reference-files.toml or delete them):
  packages/treetime/proptest-regressions/timetree/inference/__tests__/test_prop_runner.txt
```

`dev/fixtures-check` lists every tracked file under the `roots` of `dev/registry/reference-files.toml` and reports each file that no entry claims ([dev/fixtures-check#L35](../../dev/fixtures-check#L35)). The registry lists `packages/treetime/proptest-regressions` as a root ([dev/registry/reference-files.toml#L24](../../dev/registry/reference-files.toml#L24)), and claims the other three files under it with one entry each, of kind `proptest-regression` ([dev/registry/reference-files.toml#L194-L223](../../dev/registry/reference-files.toml#L194-L223)): `test_prop_resolve_polytomy.txt`, `test_prop_timetree_weak_signal.txt`, and `test_prop_marginal_dense_sparse_simulated.txt`. Each entry names its owning test, and its other fields are the same text for every proptest corpus ("auto-managed; not manually captured").

Commit `67550509e` ("date the tips that the clock filter flags from their parent") added the seed file of `packages/treetime/src/timetree/inference/__tests__/test_prop_runner.rs` without an entry. The proptest harness writes such a file by itself when a property fails, so the author commits a file that the author did not write by hand. `just check` does not run `fixtures-check`: the recipe is in the `checks_format` group ([justfile#L123](../../justfile#L123)), which only `just check-all` and its CI job run, so the routine gate of a change does not report the missing entry.

## Fix direction

Add an entry `proptest_timetree_runner` of kind `proptest-regression` for the file, with `test_prop_runner.rs` as its owning test and the same fields as the three other proptest entries.

> [!IMPORTANT]
> **Decision required.** The same failure returns with the next proptest seed file that a commit adds. Two options remove the class:
>
> - **Rule in the check**: `dev/fixtures-check` claims every `*.txt` under `proptest-regressions/` whose owning test exists at the mirrored path under `src/` (`proptest-regressions/<path>.txt` belongs to `src/<path>.rs`), and reports a seed file without such a test. The three per-file entries are then removed
> - **Keep per-file entries**: add `fixtures-check` to the fast gate of `just check`, so the author of the commit sees the missing entry
