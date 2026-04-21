# resolve_polytomies_with_options ignores merge_compressed parameter

`resolve_polytomies_with_options()` in the timetree optimization module documents `merge_compressed` as a behavior-changing parameter but never reads its value. The parameter is accepted in the function signature and passed through from callers, but the function body does not branch on it.

## Impact

The `merge_compressed` flag is intended to control whether compressed (Fitch-identical) branches are merged during polytomy resolution. Since the flag is ignored, polytomy resolution always uses the same behavior regardless of the setting, which deviates from the documented contract.

## Fix

Either implement the `merge_compressed` branching logic or remove the parameter from the function signature and its callers. If removed, document whether the current (always-on or always-off) behavior is the intended default.
