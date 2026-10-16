# Indel policy variants lack output tests

The two indel policies of branch-length optimization (`IndelPolicy::Include`, `IndelPolicy::Exclude`, selected by `--no-indels`) have unit tests of their terms in [packages/treetime/src/optimize/__tests__/test_indel_policy.rs](../../packages/treetime/src/optimize/__tests__/test_indel_policy.rs), but no test compares the optimized branch lengths of each policy against an external reference.

## Evidence

- The timetree runner golden masters are `#[ignore]`d, so they do not run
- `fn test_gm_optimize()` in [packages/treetime/src/optimize/__tests__/test_gm_optimize.rs](../../packages/treetime/src/optimize/__tests__/test_gm_optimize.rs) compares both policies with one v0 fixture. v0 has no indel term, so the fixture is a valid oracle only for the `Exclude` policy, and the `Include` policy passes only because the tolerance admits its difference
- The total branch length comparison uses `max_relative = 0.05` and the per-branch comparison uses `max_relative = 0.10`. Both tolerances predate the indel policy and are loose, so a regression of a few percent in either policy passes

> [!IMPORTANT]
> **Investigation required.** Find an independent oracle for the `Include` policy. Candidates to evaluate: a closed-form optimum for an edge whose only information is the indel count (the Poisson estimate $k/\mu$), a numerical maximization of the summed substitution and indel log-likelihood computed outside the optimizer, and a dense-versus-sparse comparison of the same policy. For `Exclude`, the existing v0 fixture applies at a tolerance that the evidence supports; measure the actual difference before tightening.
