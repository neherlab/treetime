# Dead `GTR::is_multi_site()` method

`GTR::is_multi_site()` at [packages/treetime/src/gtr/gtr.rs#L493-L495](../../packages/treetime/src/gtr/gtr.rs#L493-L495) checks `self.pi.shape().len() > 1`. Since `pi` is `Array1<f64>`, this always returns `false` by type constraint. The method is never called from any production or test code.

Two commented-out call sites exist at `gtr.rs:526` and `gtr.rs:546`, suggesting it was part of a planned multi-site GTR feature that was not completed.

## Impact

Negligible. Dead code that increases maintenance surface.

## Fix

Remove the method and the commented-out call sites. If multi-site GTR is needed later, it will use `GtrSiteSpecific` which has a different architecture.
