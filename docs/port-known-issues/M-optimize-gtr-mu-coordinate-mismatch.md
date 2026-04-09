# Optimizer and marginal propagation use incompatible branch-length scales when GTR mu ≠ 1

The branch-length optimizer and `update_marginal` disagree on the meaning of the branch-length
variable `t`. The optimizer evaluates `exp(eigvals * t)` — treating `t` as subs/site in
normalized eigenvalue units. `update_marginal` passes `t` to `expQt(t)`, which computes
`exp(eigvals * mu * t)` — treating `t` as time and relying on `mu` to convert to subs/site.
These are equivalent only when `mu = 1`. When the GTR is fitted from a timetree with year-scale
branch lengths, `mu ≈ clock_rate` (e.g., `~0.0007 subs/site/year` for SARS-CoV-2), and the two
computations diverge by a factor of roughly `1/mu ≈ 1400`. The mismatch shifts node profiles
toward the equilibrium distribution after the first optimization step resets branch lengths from
year-scale to subs-per-site scale, degrading or fully destroying the signal driving subsequent
iterations.

## Depends on

This issue co-occurs with
[M-optimize-negative-branch-length-validation.md](M-optimize-negative-branch-length-validation.md)
when `treetime optimize` is run on a timetree. That issue concerns the crash triggered by
negative timetree branch lengths; this issue concerns the wrong-results path that persists even
when negative lengths are absent or have been replaced. Both are triggered by the same class of
input (timetree-derived trees with year-scale branch lengths).

## Affected code locations

| Location | Role |
|---|---|
| [packages/treetime/src/gtr/gtr.rs#L440](../../packages/treetime/src/gtr/gtr.rs#L440) | `expQt(t)` — computes `exp(eigvals * mu * t)` using `exp_lt` |
| [packages/treetime/src/gtr/gtr.rs#L317](../../packages/treetime/src/gtr/gtr.rs#L317) | `expQt_with_rate(t, rate)` — same convention with per-site rate |
| [packages/treetime/src/gtr/gtr.rs#L330](../../packages/treetime/src/gtr/gtr.rs#L330) | `evolve()` — propagates profiles forward in time via `expQt` |
| [packages/treetime/src/gtr/gtr.rs#L383](../../packages/treetime/src/gtr/gtr.rs#L383) | `propagate_profile()` — propagates likelihoods upward via `expQt` |
| [packages/treetime/src/commands/optimize/optimize_eval.rs#L29](../../packages/treetime/src/commands/optimize/optimize_eval.rs#L29) | `evaluate_site_contributions` — computes `exp(eigvals * t)` without mu |
| [packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L21](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs#L21) | Dispatches `evaluate_site_contributions` with `gtr.eigvals` and raw `branch_length` |
| [packages/treetime/src/commands/optimize/run.rs#L152-L153](../../packages/treetime/src/commands/optimize/run.rs#L152-L153) | Main optimization loop — calls `update_marginal` then `run_optimize_mixed` in sequence |
| [packages/treetime/src/commands/optimize/optimize_unified.rs#L720](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L720) | `initial_guess_mixed` — sets BL = `subs_count / effective_length`, implicitly assuming `mu = 1` |

## Background: two incompatible conventions for `t`

`GTR::new` normalizes `W` so that `pi' * W * pi = 1` (average rate = 1 in normalized units),
and stores `mu = mu_input * average_rate_raw`. The matrix exponential is then defined as:

```
expQt(t) = V * diag(exp(eigvals * mu * t)) * V_inv
```

where `eigvals` are the eigenvalues of the *normalized* W. By convention, `t` is a physical time
parameter (e.g., years) and `mu` carries the units (subs/site/year). The product `mu * t` is
dimensionless and represents expected substitutions per site.

The optimizer uses a different convention. `evaluate_site_contributions` computes:

```
exp_ev = exp(eigvals * branch_length)
```

There is no `mu` factor. The variable named `branch_length` here means something closer to
subs/site: the optimizer's objective is self-consistent when the ML optimal `t*` is in subs/site
units. This is consistent with `expQt(t)` only if `mu = 1`.

**When mu = 1**, both callers agree: one unit of `t` equals one expected substitution per site,
and the two expressions produce the same transition probabilities for the same numeric `t`.

**When mu ≠ 1**, the same numeric value of `t` is interpreted differently:
- `expQt(t)` evaluates the matrix transition at `mu * t` subs/site
- The optimizer evaluates the matrix transition at `t` subs/site

A numeric value of `t = 0.0007` represents:
- `expQt`: `0.0007 * 0.0007 = 4.9e-7` subs/site → essentially no evolution
- optimizer: `0.0007` subs/site → correct for sc2

## Failure mechanism

Running `treetime optimize` on a SARS-CoV-2 tree produced by timetree inference:

1. **Input**: branch lengths in years, `t ≈ 1.0` for a one-year branch.

2. **GTR inference** (`infer_gtr_dense` / `infer_gtr_sparse`): accumulates `nij` (expected
   substitution counts) and `Ti` (time in state) from `update_marginal`'s stored `mut_stack`.
   `Ti` is proportional to `BL * state_freq`, so `Ti ~ 1.0 * n_sites`. The ratio
   `nij / Ti ≈ 0.0007 subs/site/year` gives `mu ≈ 0.0007` in the inferred GTR.

3. **`--branch-length-initial-guess=always`**: `initial_guess_mixed(overwrite = true)` replaces
   every edge with `t = subs_count / alignment_length ≈ 0.0007`. This is the ML subs/site
   estimate under `mu = 1`. The computation is numerically correct for that assumption but
   produces a value that is wrong for `update_marginal`'s convention.

4. **Main loop, iteration 1**:
   - `update_marginal(t = 0.0007)` calls `expQt(0.0007)`, computing
     `exp(eigvals * 0.0007 * 0.0007) = exp(eigvals * 4.9e-7)`.
   - At `4.9e-7` expected subs/site, the transition matrix is essentially identity.
   - `evolve` and `propagate_profile` return profiles that are near-indistinguishable from
     the input. No within-tree signal propagates. Leaf likelihoods wash out into uniform-ish
     distributions at every internal node.
   - `run_optimize_mixed` receives flat coefficient matrices. With
     `coefficients ≈ const`, the log-likelihood `ln(sum_c k_ic exp(lambda_c t))` is nearly
     independent of `t`, and the ML step pushes branch lengths toward 0.
   - Subsequent iterations start from near-zero branch lengths, collapsing profiles further.

5. **Result**: all branch lengths converge toward 0 or remain near-arbitrary small values.
   The output tree is numerically meaningless.

**`--branch-length-initial-guess=auto` (default, timetree input):**
`initial_guess_mixed(overwrite = false)` preserves the original year-scale BLs (`t ≈ 1.0`)
because they are finite and positive. Iteration 1 calls `update_marginal(t = 1.0)`, which
correctly computes `expQt(1.0) = exp(eigvals * 0.0007)`. Profiles are correct; the optimizer
converges properly and sets `t ≈ 0.0007`. From iteration 2 onward, `update_marginal(t = 0.0007)`
faces the same collapse as `always` mode. Whether one good iteration is enough for practical
convergence depends on dataset and `max_iter` and requires investigation (see below).

## Reproducible example

```bash
./dev/docker/run ./dev/dev r treetime -- optimize \
  --tree data/sc2/2844/tree.nwk \
  --aln data/sc2/2844/aln.fasta.xz \
  --dense false \
  --branch-length-initial-guess=always \
  --outdir tmp/p2-mismatch \
  -vvv
```

Expected (correct): branch lengths distributed in `[~5e-4, ~2e-3]` matching v0 output.
Observed (buggy): branch lengths converge toward 0 or produce wrong subs/site values.

Compare v0 reference:
```bash
./dev/docker/python treetime ancestral \
  --tree data/sc2/2844/tree.nwk \
  --aln data/sc2/2844/aln.fasta \
  --outdir tmp/p2-v0-reference
```

## v0 comparison

v0 does not expose this issue because it always works in subs/site units by the time
branch-length optimization runs. After fitting the molecular clock, v0 rescales all branch
lengths by the clock rate before ancestral reconstruction and GTR inference. As a result, the GTR
is fitted with subs-per-site branch lengths, producing `mu ≈ 1`. With `mu = 1`, `expQt(t)` and
`exp(eigvals * t)` are numerically identical for any `t`, and the coordinate mismatch does not
arise.

v1 takes the raw input tree without rescaling branch lengths to subs/site before GTR inference,
so `mu` reflects the molecular clock rate when the input is a timetree. The mismatch is a v1
regression with no v0 counterpart.

## User-facing workarounds

There is no clean workaround at the CLI level:

- **Avoid `--branch-length-initial-guess=always` on timetree input** until this is fixed. The
  `auto` default preserves the first iteration's signal and may produce acceptable results for
  well-sampled trees with `max_iter = 1`.
- **Pre-normalize the input tree**: multiply all branch lengths by the estimated clock rate
  before calling `treetime optimize`. This converts year-scale BLs to subs/site and makes GTR
  inference yield `mu ≈ 1`. Requires an external clock rate estimate (e.g., from `treetime
  clock`).
- **Use `--branch-length-initial-guess=auto` with `--max-iter=1`**: single iteration uses the
  correct year-scale profiles from `update_marginal` and produces one round of valid optimization.
  Multi-iteration results are unreliable.

## Proposed solution

The optimizer and `update_marginal` need to agree on the units of `t`. There are two approaches:

**A. Normalize GTR mu to 1 after inference.** After each GTR inference step, compute a
weighted-average scale across partitions and rescale all GTR models and all branch lengths:

```
scale = Σ_p (sequence_length_p * gtr_p.mu) / Σ_p sequence_length_p
for each partition p:   gtr_p.mu /= scale
for each edge:          branch_length *= scale
```

After renormalization, `expQt(t) = exp(eigvals * 1.0 * t) = exp(eigvals * t)`, matching the
optimizer's convention. Branch lengths are now in subs/site units. This is what v0 effectively
does by working in subs/site from the start.

**B. Include mu in the optimizer evaluator.** Multiply the eigenvalues by `gtr.mu` before
passing to `evaluate_site_contributions`, or add a `mu` parameter to the function:

```rust
let exp_ev = (eigvals * gtr.mu * branch_length).mapv(f64::exp);
```

This makes the optimizer's `t*` a time parameter (years), consistent with `update_marginal`.
Branch lengths output by the optimizer would then be in year-scale units, which changes their
interpretation relative to what users and downstream code currently expect.

Option A is closer to v0 behavior and avoids changing the semantics of the optimizer's `t`. The
rescaling step belongs after each GTR inference call, immediately before the next
`run_optimize_mixed` call. Option B is a smaller code change but has broader semantic
consequences. The right choice depends on how downstream code and output formats interpret the
branch lengths (see investigation items below).

## Investigation needed

Before implementing, the following must be confirmed:

1. **Measure actual mu values.** Run `sc2/2844`, `ebola/20`, and `tb/20` with
   `--branch-length-initial-guess=always -vvv` and capture `gtr.json` for each. Record the
   inferred `mu`. Confirm whether it is systematically far from 1 for timetree inputs across
   datasets, and whether it is close to 1 for non-timetree inputs. This grounds the severity
   assessment.

2. **`auto` mode multi-iteration behavior.** Run `--branch-length-initial-guess=auto` with
   `--max-iter 1`, `--max-iter 2`, `--max-iter 5` on `sc2/2844`. Compare final branch-length
   distributions against v0. Establish whether convergence degrades from iteration 1 to
   iteration 2 and whether the issue affects practically all timetree inputs or only `always`
   mode.

3. **Does `infer_gtr_impl` renormalize mu internally?** Read `infer_gtr_impl` in
   `packages/treetime/src/gtr/infer_gtr/common.rs` and trace through the mu computation. Check
   whether it imposes any constraint on the output mu, or whether mu is always proportional to
   the observed-rate / Ti ratio (which equals clock_rate for year-scale timetrees). The answer
   determines whether the mismatch is guaranteed to occur for any timetree input.

4. **Sparse vs. dense equivalence.** The repro uses `--dense false`. Confirm that
   `evaluate_sparse_contribution` dispatches to the same `evaluate_site_contributions` function
   and therefore omits mu identically. If so, the fix applies uniformly to both paths.

5. **Output branch-length semantics.** Determine what downstream code (Newick writer, diversity
   metrics, summary statistics) expects the branch lengths to mean after `treetime optimize`.
   If the expectation is subs/site, option A is correct. If years are expected, option B is
   correct. This governs which fix to implement.

6. **Timetree integration.** The timetree command calls its own internal optimize loop. Check
   whether the same mismatch exists there, or whether timetree's internal loop happens to work
   in a coordinate system where mu = 1 at the point of optimization.

## Tests needed

**Regression repros (currently producing wrong results with `always` mode):**
- `sc2/2844` with `--dense false --branch-length-initial-guess=always`: final branch lengths
  should be in subs/site range (~5e-4 to 2e-3), not near 0
- Same with `--dense true`: same expected range
- `ebola/20` with `--branch-length-initial-guess=always`: check branch lengths within
  expected subs/site range for ebola rates

**Consistency checks after fix:**
- For any dataset: `update_marginal` profiles after optimization should differ meaningfully from
  root priors on leaf-adjacent edges (non-collapse check — profiles at leaves must carry
  alignment signal, not be near-uniform)
- GTR mu after the renormalization step (option A) should be `1.0 ± 1e-10` for all partitions
- Branch lengths output by the optimizer should be consistent between `auto` and `always` modes
  when starting from the same initial tree

**No-regression checks:**
- Non-timetree input (e.g., a tree with BLs already in subs/site): `always` mode should still
  produce correct results after fix
- Single-iteration `auto` mode result should not change (it was already correct)
- `--branch-length-initial-guess=never` with a valid subs/site tree should not regress
