# Auspice trait reader contract is unverified after the tree I/O refactor

> [!IMPORTANT]
> **Investigation required.** The current reader has no evidence of incorrect behavior and no evidence of a verified contract. Collect the evidence listed under "Evidence required": current `AuspiceRead` implementations, malformed-value handling, and the accepted shapes from the supported Auspice schema.

The removed `tree_ir` Auspice reader silently dropped malformed trait objects and non-numeric confidence entries. The current generic reader in [`packages/treetime-io/src/auspice.rs`](../../packages/treetime-io/src/auspice.rs) delegates node conversion through `AuspiceRead`, while the repository currently contains no command-level implementation that demonstrates how typed trait values, confidence maps, and entropy are consumed.

The historical failure therefore does not establish current incorrect behavior. The remaining problem is an unverified external-format contract.

## Evidence required

- Identify every current `AuspiceRead` implementation and whether trait fields are reachable through a public command.
- Verify absent traits remain valid while present malformed values cannot disappear silently.
- Establish the accepted shapes and finite-number requirements for `value`, `confidence`, and `entropy` from the supported Auspice schema.
- Add whole-document tests only after the current conversion boundary and intended typed representation are identified.

## Readiness

Implementation waits for the evidence above. The removed `tree_ir` paths and `TreeIrTrait` type cannot be used as instructions for the current reader.
