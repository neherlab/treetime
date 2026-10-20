# Indel policy variants lack output tests

The two indel policies of branch-length optimization (`IndelPolicy::Include`, `IndelPolicy::Exclude`, selected by `--no-indels`) have unit tests of their terms in [packages/treetime/src/optimize/__tests__/test_indel_policy.rs](../../packages/treetime/src/optimize/__tests__/test_indel_policy.rs), but no test compares the optimized branch lengths of each policy against an external reference.

## Evidence

- The timetree runner golden masters are `#[ignore]`d, so they do not run
- `fn test_gm_optimize()` [packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L23-L47](../../packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L23-L47) compares the total branch length of both policies with one v0 fixture at `max_relative = 0.05`. v0 has no indel term, so the fixture is a valid oracle only for the `Exclude` policy, and the `Include` policy passes only because the tolerance admits its difference
- `fn test_gm_optimize_per_branch()` [packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L50-L90](../../packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L50-L90) compares each branch length at `max_relative = 0.10`, but it is `#[ignore]`d ([M-optimize-gm-per-branch-divergence.md](M-optimize-gm-per-branch-divergence.md)) and runs only the `Include` policy, so no per-branch comparison runs for either policy
- Both tolerances predate the indel policy and are loose, so a regression of a few percent in either policy passes

> [!IMPORTANT]
> **Investigation required.** Find an independent oracle for the `Include` policy. Candidates to evaluate: a closed-form optimum for an edge whose only information is the indel count (the Poisson estimate $k/\mu$), a numerical maximization of the summed substitution and indel log-likelihood computed outside the optimizer, and a dense-versus-sparse comparison of the same policy. For `Exclude`, the existing v0 fixture applies at a tolerance that the evidence supports; measure the actual difference before tightening.
