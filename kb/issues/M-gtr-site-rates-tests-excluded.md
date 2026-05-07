# GTR site-rate property tests excluded from build

The file [packages/treetime/src/gtr/**tests**/test_prop_gtr_site_rates.rs](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_site_rates.rs) (238 lines, 7 property tests) is not registered in [packages/treetime/src/gtr/**tests**/mod.rs](../../packages/treetime/src/gtr/__tests__/mod.rs). The `mod.rs` file lists 12 test modules but does not include `mod test_prop_gtr_site_rates;`. The file exists on disk but is silently excluded from compilation and test collection.

## Impact

The per-site rate paths in `GTR::expQt_with_rate`, `GTR::evolve`, and `GTR::propagate_profile` (at [packages/treetime/src/gtr/gtr.rs#L319](../../packages/treetime/src/gtr/gtr.rs#L319), [gtr.rs#L356-L367](../../packages/treetime/src/gtr/gtr.rs#L356-L367), [gtr.rs#L405-L415](../../packages/treetime/src/gtr/gtr.rs#L405-L415)) are wired to production via [packages/treetime/src/representation/partition/marginal_helpers.rs](../../packages/treetime/src/representation/partition/marginal_helpers.rs) and have no test coverage in the current build. The only existing coverage for per-site rate behavior is through `test_prop_gtr_site_specific.rs`, which tests the separate `GTRSiteSpecific` class, not the scalar `GTR` per-site branch.

## Affected code

- Missing registration: [packages/treetime/src/gtr/**tests**/mod.rs](../../packages/treetime/src/gtr/__tests__/mod.rs) -- `test_prop_gtr_site_rates` not listed
- Excluded tests: [packages/treetime/src/gtr/**tests**/test_prop_gtr_site_rates.rs](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_site_rates.rs)

## Fix

Add `mod test_prop_gtr_site_rates;` to [packages/treetime/src/gtr/**tests**/mod.rs](../../packages/treetime/src/gtr/__tests__/mod.rs). Verify the tests compile and run. If any tests fail, investigate and file the failure as a separate known issue rather than leaving the tests excluded.
