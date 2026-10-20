# Timetree with input branch lengths aborts on a point-by-point division

`treetime timetree --branch-length-mode=input` fails on every dataset and writes no output:

```
Error:
   0: Cannot divide point by point: operation not well-defined

Location: packages/treetime-distribution/src/distribution_ops/divide.rs:34
```

## Symptom and reproduction

```bash
./dev/docker/run just r treetime timetree --branch-length-mode=input \
  --tree=data/flu/h3n2/20/tree.nwk \
  --metadata=data/flu/h3n2/20/metadata.tsv \
  --alignment=data/flu/h3n2/20/aln.fasta.xz \
  --output-all=tmp/repro-input-bl
```

The same abort occurs on `data/ebola/20`, `data/zika/20`, and on `data/ebola/362` with `--keep-root`. These are the `--branch-length-mode=input` cases of `dev/smoke` (the `branch-input`, `branch-input-all-outputs` and `input-keeproot` variants of `timetree`), which the matrix declares as expected failures [dev/smoke.toml#L659-L679](../../dev/smoke.toml#L659-L679). The same abort occurs on the simulated data of [validation/README.md](../../validation/README.md). The test `test_pipeline_input_mode_undated_single_child_root_fails_with_point_division_known_issue_h_timetree_input_branch_lengths_abort_on_point_division` asserts the failure [packages/treetime/src/timetree/__tests__/test_pipeline.rs#L212](../../packages/treetime/src/timetree/__tests__/test_pipeline.rs#L212).

Input mode reconstructs no ancestral sequences, so `--output-all` skips the reconstructed FASTA, and `--output-reconstructed-nuc-fasta` fails before inference.

## Impact and scope

- `--branch-length-mode=input` is unusable for `timetree`: the command produces no output at all. The default mode (`marginal`) completes on the same datasets
- These cases are excluded from every before/after output comparison, so this part of the timetree pipeline has no regression coverage
- Before the forward-pass division existed, the same cause showed as a different symptom: the Nexus output carried dates only on the tips (19 of 37 nodes on `data/flu/h3n2/20`), because the product of two points at slightly different times is empty

## Mechanism

`fn create_branch_distributions_input_mode()` builds every branch distribution as a degenerate `Distribution::point(time_length, 0.0)` [packages/treetime/src/timetree/inference/runner.rs#L308-L332](../../packages/treetime/src/timetree/inference/runner.rs#L308-L332). A `TODO` there names the v0 interpolator as the missing piece.

Leaf time distributions are points at the sampling dates. The backward pass convolves each child point with the negated branch point and gets a point message for the parent. The forward pass then forms the cavity distribution: it divides the parent posterior by the message from the child [packages/treetime/src/timetree/inference/forward_pass.rs#L192-L197](../../packages/treetime/src/timetree/inference/forward_pass.rs#L192-L197). `fn distribution_division()` defines division only by a `Function` and rejects a `Point` or `Range` divisor [packages/treetime-distribution/src/distribution_ops/divide.rs#L31-L34](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L31-L34), so the pass aborts.

A point product has the same weakness: `fn multiply_point_point()` returns `Empty` when the two times differ by more than `1e-9` [packages/treetime-distribution/src/distribution_ops/multiply.rs#L53-L62](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L53-L62), and floating-point accumulation through the tree gives differences of `1e-7` to `1e-5`.

### Why v0 does not fail

v0's `BranchLenInterpolator` never uses delta distributions for branch lengths, also in input mode. It builds a distribution over a grid [packages/legacy/treetime/treetime/branch_len_interpolator.py#L64-L102](../../packages/legacy/treetime/treetime/branch_len_interpolator.py#L64-L102):

- branch length below 0.05: Poisson likelihood of the number of mutations
- longer branches: Gaussian, with a variance corrected for saturation

v0's `Distribution.multiply()` evaluates the other factors at the position of a single delta and raises `ArithmeticError` for more than one delta [packages/legacy/treetime/treetime/distribution.py#L82-L149](../../packages/legacy/treetime/treetime/distribution.py#L82-L149). Because v0 never puts two deltas into one product, it never multiplies two points at different times.

## Scientific context

Branch lengths are evolutionary distances in expected substitutions per site. A molecular clock converts them to time: `time = branch_length / clock_rate`.

**Poisson for short branches (below 0.05 substitutions per site).** Substitutions occur as a Poisson process along a branch. For a branch of true length `k` on a sequence of length `L`, the expected number of substitutions is `k*L`, and the probability of `n` substitutions is

```
P(n | k, L) = exp(-k*L) * (k*L)^n / n!
```

Given the observed branch length `l`, the likelihood of a candidate time span `dt`, which implies `k = dt * clock_rate`, is

```
log P(k | l, L) = -k*L + l*L*log(k)
```

**Gaussian for long branches (0.05 substitutions per site and above).** On longer branches a site can mutate more than once, so the observed difference `p` saturates at `p0 = 1 - sum(Pi^2)` (about 0.75 for nucleotides with uniform frequencies). The observed difference relates to the true distance `l` as

```
p = p0 * (1 - exp(-l/p0))
```

Finite sampling gives `p` the variance `p*(1-p)/L`. Through the saturation correction, the variance of the branch length estimate becomes

```
sigma^2 = p0 * (exp(l/p0) - 1) * (exp(l/p0) - p0*(exp(l/p0) - 1)) / L
```

A Gaussian centered at the observed branch length with this variance captures the uncertainty from finite sampling and saturation.

**Why points fail.** A point distribution has no uncertainty. When the backward and forward passes arrive at slightly different times, the product of two points is empty and the quotient is undefined, which breaks the message passing.

## Fix

Reimplement `create_branch_distributions_input_mode()` as v0's `BranchLenInterpolator` input mode:

- Poisson (v0 lines 78-86): `log_prob = -k*L + l*L*log(k)`, with `k` the candidate branch length, `l` the observed branch length, and `L` the sequence length
- Gaussian (v0 lines 87-102): `sigma_sq = p0*(exp(l/p0)-1)*(exp(l/p0)-p0*(exp(l/p0)-1))*one_mutation`, with `p0 = 1 - sum(Pi^2)`

v1 has no production Poisson branch-length distribution. A golden-master test helper builds one: `create_poisson_branch_distributions()` [packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_support.rs#L125](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_support.rs#L125). The sequence length is available: without an alignment, v1 already requires `--sequence-length`, as v0 does [packages/treetime/src/timetree/pipeline.rs#L236-L248](../../packages/treetime/src/timetree/pipeline.rs#L236-L248).

## Validation

- The `--branch-length-mode=input` cases of `dev/smoke` and of the validation suite complete, and their expected failures are removed from the matrix
- The per-edge distribution in input mode matches v0's on a small tree
- Every internal node of `data/flu/h3n2/20` gets a date in the Nexus output
- The test that asserts the abort is replaced by one that asserts the dates

## Related issues

- [M-timetree-input-mode-edges-without-branch-length-carry-no-evidence.md](M-timetree-input-mode-edges-without-branch-length-carry-no-evidence.md)
- [N-timetree-gtr-json-missing-input-bl.md](N-timetree-gtr-json-missing-input-bl.md)
- [M-timetree-weak-signal-errors-differ-by-code-path.md](M-timetree-weak-signal-errors-differ-by-code-path.md)
