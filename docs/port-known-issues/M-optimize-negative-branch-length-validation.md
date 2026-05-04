# Negative branch lengths not validated or reported across optimizer modes

The branch length optimizer does not validate, warn about, or uniformly handle negative input branch lengths. Depending on which `--branch-length-initial-guess` mode is used and whether the affected edges carry indels, the consequences range from silent wrong results to a debug-build panic, with no message to the user in any case.

## Background

Timetree inference (including TreeTime's own marginal timetree pipeline) can produce trees where a small number of edges have slightly negative branch lengths (e.g. `-0.000189`). This artifact occurs because marginal tip-date inference optimizes each node time independently; rounding or numerical constraint relaxation can place a child node infinitesimally earlier than its parent. These values are meaningless for substitution-rate optimization: the Poisson indel log-likelihood $\ell(t) = k \ln(\mu t) - \mu t$ is defined only for $t > 0$ when $k > 0$, and substitution matrix exponentiation $\exp(Q t)$ with $t < 0$ is mathematically valid but physically nonsensical.

## Affected code locations

| Location                                                                                                                                             | Role                                                                       |
| ---------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------- |
| [packages/treetime/src/commands/optimize/run.rs#L710-L738](../../packages/treetime/src/commands/optimize/run.rs#L710-L738)                           | Mode dispatch -- selects what to do before optimization                    |
| [packages/treetime/src/commands/optimize/run.rs#L650-L665](../../packages/treetime/src/commands/optimize/run.rs#L650-L665)                           | `is_branch_length_missing` and `any_edge_missing_branch_length` predicates |
| [packages/treetime/src/commands/optimize/optimize_unified.rs#L583-L597](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L583-L597) | Per-edge bootstrap inside `run_optimize_mixed`                             |
| [packages/treetime/src/commands/optimize/optimize_unified.rs#L696-L704](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L696-L704) | Skip condition inside `initial_guess_mixed`                                |
| [packages/treetime/src/commands/optimize/optimize_indel.rs#L29-L31](../../packages/treetime/src/commands/optimize/optimize_indel.rs#L29-L31)         | `debug_assert!(t > 0.0)` in `poisson_indel_log_lh`                         |

## Behaviour per mode

### `--branch-length-initial-guess=auto` (default)

`initial_guess_mixed` is called with `overwrite_valid = false`. The skip condition at [optimize_unified.rs#L702](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L702) is:

```rust
if bl.is_finite() && (bl > 0.0 || indel_count == 0) {
    continue;
}
```

For an edge with a **negative BL and indels present**: `bl > 0.0 = false`, `indel_count == 0 = false` → condition is `false` → edge is **rewritten** from substitution counts. This specific sub-case does not crash.

For an edge with a **negative BL and no indels**: `bl > 0.0 = false`, `indel_count == 0 = true` → condition is `true` → edge is **silently skipped**, preserving the negative value. Newton optimization then starts from a negative branch length. The Newton step is clamped so that the result is non-negative after the first iteration, but the initial gradient is evaluated at a physically nonsensical `t < 0`. No message is emitted.

In both cases: **no warning**, no indication that the input tree contained negative branch lengths.

### `--branch-length-initial-guess=never`

`initial_guess_mixed` is not called. The only guard is:

```rust
fn is_branch_length_missing(bl: Option<f64>) -> bool {
    bl.is_none_or(|v| v.is_nan())
}
```

This accepts any finite value, including negatives. The error message in `Never` mode says "fails if any edge has a missing or invalid branch length" (from the `--help` text), but **negative values are not considered invalid** by this predicate.

Negative branch lengths pass silently into `run_optimize_mixed`. The per-edge bootstrap at [optimize_unified.rs#L583](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L583) handles `branch_length == 0.0` only, not `< 0.0`:

```rust
if branch_length == 0.0 && indel_count > 0 { … }
```

A negative value reaches `evaluate_with_indels → poisson_indel_log_lh` with `t < 0` and `k > 0`:

- **Debug build**: `debug_assert!(t > 0.0, "poisson_indel_log_lh requires t > 0 when k > 0, got t={t}")` - **panic/crash**.
- **Release build**: The assert is compiled out. `ln(mu * t)` with `t < 0` computes `ln(negative)` = `NaN`, which propagates silently through the rest of optimization. Results are silently wrong.

For edges with no indels in `Never` mode, the substitution-side evaluation (`exp(eigvals * t)` with `t < 0`) produces mathematically valid but physically wrong exponentials. No crash, but the optimizer starts from an incorrect position. **No message emitted.**

### `--branch-length-initial-guess=always`

`initial_guess_mixed` is called with `overwrite_valid = true`, which unconditionally rewrites every edge from substitution counts regardless of the current value. Negative BLs are replaced by non-negative substitution-rate estimates. The crash and numerical problem do not occur.

**No warning is emitted.** The user has no feedback that their input contained negative branch lengths.

## Summary table

| Mode              | Negative BL (+ indels)                | Negative BL (no indels)                      | User notified? |
| ----------------- | ------------------------------------- | -------------------------------------------- | -------------- |
| `auto`            | Rewritten silently (no crash)         | Preserved silently, Newton starts at `t < 0` | No             |
| `always`          | Rewritten silently (no crash)         | Rewritten silently (no crash)                | No             |
| `never` (debug)   | **Panic**                             | Wrong starting point                         | No             |
| `never` (release) | NaN propagation, silent wrong results | Wrong starting point                         | No             |

## Inconsistencies

Four separate validity predicates exist across the codebase. They disagree on what constitutes an invalid branch length:

| Predicate                                           | Rejects `None` | Rejects `NaN` | Rejects `< 0`             |
| --------------------------------------------------- | -------------- | ------------- | ------------------------- |
| `is_branch_length_missing` (run.rs)                 | ✓              | ✓             | ✗                         |
| `initial_guess_mixed` skip (Auto mode, no indels)   | ✓              | ✓             | ✗ (−BL treated as valid)  |
| `initial_guess_mixed` skip (Auto mode, with indels) | ✓              | ✓             | ✓ (−BL triggers rewrite)  |
| `run_optimize_mixed` bootstrap                      | n/a            | n/a           | ✗ (handles `== 0.0` only) |

The `Never`-mode guard and the initial-guess skip condition for no-indel edges both fail to treat negative branch lengths as invalid, while the skip condition for indel-bearing edges does the right thing.

## v0 comparison

v0 does not model indels in branch length optimization. Negative branch lengths in the input tree are passed directly to the substitution log-likelihood without assertion (`log(exp(Q * t))` is evaluated, which is finite for `t < 0` under most GTR models). v0 does not warn about or reject negative branch lengths.

## User-facing workarounds

If you are running `treetime optimize` on a tree produced by a previous timetree inference and observe a panic or unexpected results:

- **`--branch-length-initial-guess=always`** (recommended for timetree input): recomputes all branch lengths from substitution counts before optimization. Eliminates the crash and the wrong-starting-point problem at the cost of discarding timetree-calibrated branch lengths.
- **`--branch-length-initial-guess=auto`** (default): safe when the tree has no indels on negative-BL edges. Still silently preserves negative values on no-indel edges and produces suboptimal (but usually close) results.
- **Avoid `--branch-length-initial-guess=never`** with timetree input: negative branch lengths are not caught, leading to a debug panic or silent NaN in release.
- **As a pre-processing step**, filter or zero-clamp negative branch lengths in the input Newick before passing to `treetime optimize`.

## Proposed solution

The core fix is one canonical validity predicate that treats negative branch lengths as invalid, used consistently across all four sites, plus a user-facing warning that fires unconditionally before the mode dispatch.

### Warning message

Before the mode dispatch, scan all edges for invalid branch lengths (`None`, `NaN`, negative) and emit a `warn!` listing each affected edge by name. Edge names can be derived from the graph: source and target node keys are obtainable from an edge reference via `.source()` / `.target()`, and node names are stored in the node payload. See how edge labels are constructed in the ancestral tests for a working pattern. The warning should suggest `--branch-length-initial-guess=always` as the path-of-least-resistance fix for users with timetree input.

### Validity predicate

Define a single canonical predicate that treats a branch length as valid only when it is present, finite, and non-negative. Replace all four ad-hoc checks with this predicate. The intent is: any value that cannot be used as a physical branch length in a likelihood expression is invalid and must be flagged.

### `Auto` skip condition in `initial_guess_mixed`

The current skip logic has two branches of intent: skip if the value is already good, rewrite otherwise. The indel-count term in the condition exists for a specific reason (zero BL is invalid when indels are present, even though zero is finite and non-negative - the Poisson derivative diverges there). Any replacement must preserve that nuance, not flatten it into a single validity check. Be careful about collapsing or simplifying this condition without reading the surrounding comment and understanding why zero is treated differently.

### `Never`-mode guard

Extend to reject negative branch lengths with the same error path, and update the error message to state what was found (mirroring the warning format above). The contract advertised in `--help` says it fails on "missing or invalid" branch lengths; negative values should match that stated contract.

### Bootstrap in `run_optimize_mixed`

The per-edge bootstrap that re-seeds from a positive starting point currently only triggers on exactly-zero branch lengths. It should also catch values below zero, as they are equally unusable in this context. This is a defensive backstop; it should not be the primary fix.

### `poisson_indel_log_lh` assert

Promoting the `debug_assert!` to a proper `Result`-returning error requires changing the function signature. That change propagates upward through `evaluate_with_indels` and `evaluate_with_indels_log_lh_only`, and from there into the Brent adapter structs and every call site inside `run_optimize_mixed`. Read the full call chain before touching this function. The change is correct and makes release builds safe, but the blast radius is larger than the single-line assert suggests.

## Tests needed

The fix touches four distinct code paths across three different entry modes. Coverage should be exhaustive across the matrix. Start with red tests that reproduce the reported failures, then expand to the adjacent cases.

**Regression repros (currently failing or crashing):**

- `Never` mode + negative BL + indels on the affected edge: should return an error, not panic
- `Never` mode + negative BL + no indels: should return an error (negative is invalid by contract)
- `Auto` mode + negative BL + no indels: should emit a warning and not preserve the negative value into optimization

**Correctness after fix:**

- `Auto` mode + negative BL + no indels: branch length is replaced, result is finite and ≥ 0
- `Auto` mode + negative BL + indels: same
- `Always` mode + any invalid BL: all replaced, warning still emitted
- `Never` mode + all valid (positive) BLs: no error, no regression
- Warning message lists the correct edges with their values
- `Never`-mode error message lists the correct edges

**Boundary conditions:**

- Exactly-zero BL on an indel-free edge (should remain valid - zero branch length is a meaningful optimizer result)
- Exactly-zero BL on an indel-bearing edge (handled by existing bootstrap, should not regress)
- Mixed tree: some negative, some zero, some positive BLs across different edge types
- Tree with no invalid BLs: no warning emitted, no regression in any mode

**`poisson_indel_log_lh` promotion (if done):**

- Calling with `t < 0` and `k > 0` returns an error in both debug and release builds
- Calling with `t = 0` and `k = 0` still returns the correct zero-indel result
