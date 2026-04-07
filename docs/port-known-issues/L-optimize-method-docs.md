# Per-edge optimization: documentation and ledger updates for 6-method selection

## Problem

The ledger documents describe the old 3-method design. Several tradeoff claims in the intentional-changes doc are incorrect (e.g. "Lost: sqrt(t) reparameterization" -- no longer lost). The algorithm, feature, and test inventories do not reflect the new methods. The convergence behavior of each method is undocumented in the code.

## Changes

### 1. Intentional changes doc

`docs/port-intentional-changes/optimize-newton-raphson-per-edge.md`

The title says "uses Newton-Raphson instead of Brent" but the default is now Brent. Rewrite to describe the 6-method selection feature. Key tradeoff updates:

- "Lost: sqrt(t) reparameterization" is now recovered -- available as `brent-sqrt` (default) and `newton-sqrt`
- New: $\ln(t)$ reparameterization (`brent-log`, `newton-log`) -- eliminates indel Hessian singularity entirely. Not present in v0 or other phylogenetic tools for branch length optimization. Precedented by coalescent Tc optimization in this codebase.
- New: method selection via `--opt-method` (6 methods)
- Still lost: regularization penalty (`exp(t^4/10000)` for marginal profiles in v0)
- Still lost: Hamming distance fallback (v0 returns Hamming distance on Brent failure; v1 returns current branch length)

### 2. Proposals doc

`docs/port-proposals/optimize-convergence-and-method-choice.md`

Update P5 status to fully implemented with 6 methods. Update the per-edge optimization checklist.

### 3. Known issues index

`docs/port-known-issues/_index.md`

Verify `M-optimize-branch-method-selection.md` row is gone (already deleted by earlier commit on this branch). Add rows for the method selection issues or remove them as they are resolved.

### 4. Feature inventory

`docs/port-feature-inventory/_index.md`

Updates (line numbers approximate, verify before editing):

- Line ~425: `[ ] sqrt(t) reparameterization` change to `[x]` -- now available as `brent-sqrt` and `newton-sqrt`
- Line ~431: `[ ] Damping in marginal loop` change to `[x]` -- implemented via `--damping` flag
- Add new entries:
  - `[x] Per-edge optimization method selection (--opt-method: 6 methods)`
  - `[x] ln(t) reparameterization (brent-log, newton-log)`
  - `[x] Brent's method for per-edge branch length optimization (brent, brent-sqrt, brent-log)`

### 5. Algorithm inventory

`docs/port-algo-inventory/optimization.md`

Add entries for the 4 new optimization methods alongside existing Newton-t and Brent-t entries. Each entry should include:

**Newton-sqrt** (`newton_sqrt_inner` in `optimize_newton.rs`):

- Algorithm: Newton-Raphson in $\sqrt{t}$ space
- Chain rule: $d\ell/ds = 2s \cdot d\ell/dt$, $d^2\ell/ds^2 = 4s^2 \cdot d^2\ell/dt^2 + 2 \cdot d\ell/dt$
- Effect: reduces indel Hessian singularity from $O(1/t^2)$ to $O(1/t)$
- Convergence: quadratic near optimum, step-size criterion
- Precedent: v0 uses $\sqrt{t}$ reparameterization with Brent

**Newton-log** (`newton_log_inner` in `optimize_newton.rs`):

- Algorithm: Newton-Raphson in $\ln(t)$ space
- Chain rule: $d\ell/du = t \cdot d\ell/dt$, $d^2\ell/du^2 = t^2 \cdot d^2\ell/dt^2 + t \cdot d\ell/dt$
- Effect: eliminates indel Hessian singularity entirely ($\ell''_{\text{indel}} = -\mu t$, bounded)
- Convergence: quadratic near optimum, step-size criterion with per-space tolerance
- Novel: not present in v0 or other phylogenetic tools for branch length optimization

**Brent-sqrt** (`brent_sqrt_inner` in `optimize_brent.rs`):

- Algorithm: Brent's method in $\sqrt{t}$ space via `argmin::BrentOpt`
- Bracket: $[\sqrt{\text{min\_bl}},\; \sqrt{\text{upper}}]$, cost function evaluates $-\ell(s^2)$
- Effect: smooths objective for better parabolic interpolation
- Convergence: order ~1.325 (parabolic interpolation), derivative-free
- Precedent: matches v0 exactly (same algorithm, same parameterization). Default method.

**Brent-log** (`brent_log_inner` in `optimize_brent.rs`):

- Algorithm: Brent's method in $\ln(t)$ space via `argmin::BrentOpt`
- Bracket: $[\ln(\text{min\_bl}),\; \ln(\text{upper})]$, cost function evaluates $-\ell(e^u)$
- Effect: smoothest objective surface of all Brent variants
- Convergence: order ~1.325, derivative-free
- Precedent: coalescent Tc optimizer uses `BrentOpt` in log-space at `packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs:62`

Distinct from existing Brent entries for clock (root optimization), coalescent Tc, and polytomy resolution -- those are different optimization targets using the same `BrentOpt` solver.

### 6. Test inventory

`docs/port-test-inventory/`

Add entry for the rewritten `test_optimize_method.rs` listing all test functions and the 5 verification criteria:

- C1: local optimality (lh at optimum exceeds neighbors at 3 delta scales)
- C2: stationarity (Newton-specific: implied step below tolerance at optimum)
- C3: cross-method log-likelihood agreement (all 6 methods agree within tolerance)
- C4: bracket validity (Brent-specific: optimum lh exceeds bracket endpoints)
- C5: cross-conditioning (newton-log >= newton-sqrt >= newton-t on indel-bearing edges)

Also note: 10 test files parameterized across all 6 methods, golden master switched to `brent-sqrt`, smoke tests added for all non-default methods.

### 7. Convergence tolerance comments

Add inline doc comments in each inner-loop function documenting what the tolerance means in both the parameterized space and $t$-space:

- `newton_inner`: step-size in $t$-space, relative 0.1% with absolute floor
- `newton_sqrt_inner`: step-size in $s$-space, maps to tighter $t$-tolerance near zero ($ds/dt = 1/(2\sqrt{t})$ amplifies precision)
- `newton_log_inner`: step-size in $u$-space, natural relative tolerance in $t$-space ($dt/t \approx du$)
- `brent_inner`: BrentOpt default tolerance in $t$-space
- `brent_sqrt_inner`: BrentOpt default tolerance in $s$-space (same tightening-near-zero property)
- `brent_log_inner`: BrentOpt default tolerance in $u$-space (same relative-tolerance property)

## Verification

No code behavior change. Verify documents render correctly and cross-references resolve.

## Dependencies

- Depends on: nothing
- Depended on by: nothing (terminal)

## Cross-references

- Intentional changes doc: `docs/port-intentional-changes/optimize-newton-raphson-per-edge.md`
- Proposals doc: `docs/port-proposals/optimize-convergence-and-method-choice.md`
- Known issues index: `docs/port-known-issues/_index.md`
- Feature inventory: `docs/port-feature-inventory/_index.md`
- Algorithm inventory: `docs/port-algo-inventory/optimization.md`
- Test inventory: `docs/port-test-inventory/`
