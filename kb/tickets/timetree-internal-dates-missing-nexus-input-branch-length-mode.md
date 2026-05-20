# Internal node dates missing in nexus for input branch length mode

`timetree.nexus` contains date annotations only on tips when using `--branch-length-mode=input`. Default mode correctly annotates all nodes (tips + internal).

## Scope

flu/h3n2/20 default: 37/37 annotations (19 tips + 18 internal). flu/h3n2/20 with `--branch-length-mode=input`: 19/37 annotations (tips only).

This affects ALL datasets when using input branch length mode. Coalescent and skyline modes follow separate code paths and are not part of this issue.

## Root cause

Point distribution multiplication fails when time coordinates differ by more than 1e-9, returning `Distribution::Empty`.

### Code path

1. `create_branch_distributions_input_mode()` at
   [`runner.rs#L143-L166`](../../packages/treetime/src/timetree/inference/runner.rs#L143) creates `Distribution::point(time_duration, 1.0)` for each edge.

2. Leaf time distributions are Point distributions at collection dates.

3. Backward pass at
   [`backward_pass.rs#L62-L77`](../../packages/treetime/src/timetree/inference/backward_pass.rs#L62):
   - Convolves child Point with negated branch Point, yielding parent Point
   - Multiplies child messages via `distribution_multiplication()`
   - Point x Point works when `|t_a - t_b| <= 1e-9`

4. Forward pass at
   [`forward_pass.rs#L69-L82`](../../packages/treetime/src/timetree/inference/forward_pass.rs#L69):
   - `dist_from_parent = convolution(parent_time_dist, branch_dist)` - Point
   - `combined = multiplication(dist_from_parent, subtree_dist)` - Point x Point
   - When these Points have slightly different `t` values (floating-point
     accumulation through tree depth), `multiply_point_point()` returns Empty

5. `multiply_point_point()` at
   [`multiply.rs#L55-L65`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L55):

   ```rust
   const EPS: f64 = 1e-9;
   if (a.t() - b.t()).abs() > EPS {
     return Ok(Distribution::empty());
   }
   ```

6. Empty distribution yields `likely_time() = None`, so internal nodes have no
   time value and no date annotation in nexus output.

### Why v0 does not have this issue

v0's `BranchLenInterpolator` at [`branch_len_interpolator.py#L64-L102`](../../packages/legacy/treetime/treetime/branch_len_interpolator.py#L64) NEVER uses delta (point) distributions for branch lengths, even in input mode. Instead, it creates a rich probability distribution over a grid:

- Short branches (< 0.05): Poisson-like distribution
- Longer branches: Gaussian approximation with variance derived from sequence
  length and branch length

v0's `Distribution.multiply()` at [`distribution.py#L82-L149`](../../packages/legacy/treetime/treetime/distribution.py#L82) handles delta functions specially:

- Exactly one delta: evaluates all other distributions at the delta's position
- More than one delta: raises `ArithmeticError`

Since v0 never creates multiple delta distributions in the multiplication chain, it never hits the case where two deltas at different positions must be multiplied.

### Scientific context

Branch lengths in phylogenetics represent evolutionary distance measured in expected substitutions per site. Converting to time requires a molecular clock model: `time = branch_length / clock_rate`.

**Why Poisson for short branches (< 0.05 subs/site)**

Substitutions occur as a Poisson process along branches. For a branch of true length `k` (subs/site) on a sequence of length `L`, the expected number of substitutions is `k*L`, and the probability of observing `n` substitutions is:

```
P(n | k, L) = exp(-k*L) * (k*L)^n / n!
```

Given observed branch length `l` (from tree inference), the likelihood of candidate time-duration `dt` (implying branch length `k = dt * clock_rate`) is:

```
log P(k | l, L) = -k*L + l*L*log(k)
```

This is the Poisson model used for short branches where each site mutates at most once.

**Why Gaussian for long branches (>= 0.05 subs/site)**

For longer branches, the same site can mutate multiple times (saturation). Observed differences `p` between sequences saturate at `p0 = 1 - sum(Pi^2)` (~0.75 for nucleotides with uniform frequencies). The relationship between observed `p` and true distance `l` is:

```
p = p0 * (1 - exp(-l/p0))
```

The variance of observed `p` due to finite sampling is `p*(1-p)/L`. Propagating through the non-linear saturation correction, the variance in branch length estimate becomes:

```
sigma^2 = p0 * (exp(l/p0) - 1) * (exp(l/p0) - p0*(exp(l/p0) - 1)) / L
```

A Gaussian centered at the observed branch length with this variance captures the uncertainty from both finite sampling and saturation correction.

**Why Point distributions fail**

Point (delta) distributions assume zero uncertainty. When backward and forward passes compute slightly different time estimates (floating-point accumulation), multiplying two Points at different locations returns Empty. This breaks the message-passing algorithm.

## Repro

```bash
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 \
  --branch-length-mode=input \
  --tree=data/flu/h3n2/20/tree.nwk \
  --dates=data/flu/h3n2/20/metadata.tsv \
  --outdir=tmp/repro-input-bl data/flu/h3n2/20/aln.fasta.xz
grep -oP '\[&[^\]]*\]' tmp/repro-input-bl/timetree.nexus | wc -l
# 19 (expected: 37)
```

## Diagnostic

Add logging to `multiply_point_point()` to observe the epsilon failures:

```rust
if (a.t() - b.t()).abs() > EPS {
  log::debug!(
    "multiply_point_point: t_a={:.12e} t_b={:.12e} diff={:.12e}",
    a.t(), b.t(), (a.t() - b.t()).abs()
  );
  return Ok(Distribution::empty());
}
```

Typical output shows differences in the 1e-7 to 1e-5 range, well above the 1e-9 threshold but small enough to be floating-point accumulation rather than actual time disagreement.

## Fix

Reimplement `create_branch_distributions_input_mode()` to match v0's `BranchLenInterpolator` input mode ([`branch_len_interpolator.py#L64-L102`](../../packages/legacy/treetime/treetime/branch_len_interpolator.py#L64)):

```
if branch_length < 0.05:
    Poisson distribution
else:
    Gaussian approximation with saturation correction
```

v1 Poisson already implemented: [`utils.rs#L77-L118`](../../packages/treetime/src/timetree/utils.rs#L77) `create_poisson_branch_distributions()` - reuse or inline.

v0 reference implementations:

- Poisson (lines 78-86): `log_prob = -kL + lL*log(k)` where k is candidate
  branch length, l is observed branch length, L is sequence length
- Gaussian (lines 87-102): variance `sigma_sq = p0*(exp(l/p0)-1)*(exp(l/p0)-p0*(exp(l/p0)-1))*one_mutation`
  where `p0 = 1 - sum(Pi^2)` accounts for saturation

Sequence length: v0 requires either `--aln` or `--sequence-length` ([`wrappers.py#L381-383`](../../packages/legacy/treetime/treetime/wrappers.py#L381)). v1 should enforce the same constraint.

## Related issues

- Source: [M-timetree-internal-dates-missing-input-bl.md](../issues/M-timetree-internal-dates-missing-input-bl.md) -- delete after full resolution

This input-branch-length variant has a distinct root cause from the other timetree missing-date failures. Here the problem is Point-distribution multiplication tolerance, not branch-grid resolution or normalized product handling.

That distinction matters because the symptom is similar - internal nodes lose dates - but the remedy is different. Fixes in the grid-based branches do not address the Point x Point failure in input branch length mode.
