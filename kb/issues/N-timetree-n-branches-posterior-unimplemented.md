# Point-time lineage counts of the coalescent prior are unimplemented

`timetree` has no `--n-branches-posterior` flag until a choice of lineage counts exists. In v0 the flag is a boolean (`action='store_true'`, [packages/legacy/treetime/treetime/argument_parser.py#L311-L318](../../packages/legacy/treetime/treetime/argument_parser.py#L311-L318)). It selects how `Coalescent.calc_branch_count()` computes the lineage count $k(t)$: from the point times of the nodes (default), or from their marginal posterior time distributions [packages/legacy/treetime/treetime/merger_models.py#L84-L114](../../packages/legacy/treetime/treetime/merger_models.py#L84-L114).

The v1 coalescent prior always computes expected lineage counts from node time distributions: `fn LineageExpectation::from_gaussian()` from the Gaussian node marginals in the init phase [packages/treetime/src/coalescent/mean_field.rs#L309](../../packages/treetime/src/coalescent/mean_field.rs#L309), and `fn LineageExpectation::from_posterior()` from the survival function of each node's posterior afterwards, with a point time only for a node without a posterior distribution [packages/treetime/src/coalescent/mean_field.rs#L296-L307](../../packages/treetime/src/coalescent/mean_field.rs#L296-L307). The v0 point-time default has no v1 counterpart.

> [!IMPORTANT]
> **Decision required.** Whether to add the v0 point-time lineage counts as a choice, with the boolean flag `--n-branches-posterior` that selects the posterior counts of v1. Evidence: the v1 prior has no point-time lineage counts (see [kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md)), and v0 uses them by default.

The CLI rejects the v0 spelling with the message for unavailable flags: the entry `--n-branches-posterior` with `Then::Unavailable` in [packages/app-cli/src/cli/v0_args.rs](../../packages/app-cli/src/cli/v0_args.rs).

## Proposed implementation plan

The plan stops at the open decision above. One short plan per option:

- **Option A: no point-time counts.** Keep the v1 posterior counts as the only choice. Steps: change the `--n-branches-posterior` entry in `packages/app-cli/src/cli/v0_args.rs` from `Then::Unavailable` to a message that the v1 prior always uses posterior lineage counts, so a v0 script that passes the flag gets the same meaning; record the v0 point-time default as a difference from v0 in `kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md` (needs user approval). Tests: one case in `packages/app-cli/src/cli/__tests__/test_v0_args.rs` that asserts the new message with `assert_err!`. Then delete this issue
- **Option B: add point-time counts as a choice.** Steps:
  1. Add a lineage-count enum (point, posterior) to the coalescent prior options of `timetree`, with a clap flag and a config key in kebab-case; choose the v1 default with the user
  2. Add `fn LineageExpectation::from_point_times()` in `packages/treetime/src/coalescent/mean_field.rs` that builds `NodeSurvival::point(time)` for every node with a time, through `fn from_survivals()`, and select it in place of `from_posterior()` and `from_gaussian()` when the point mode is set
  3. Map the v0 flag `--n-branches-posterior` in `v0_args.rs` to the posterior mode, and v0 calls without the flag to the point mode only if the user approves that v0 default
  4. Tests: a unit test of `from_point_times()` on a small hand-made tree whose step function $k(t)$ is counted by hand (oracle: the lineage count definition of `Coalescent.calc_branch_count()` in `packages/legacy/treetime/treetime/merger_models.py#L84-L114` with point times); a v0 oracle comparison of the coalescent log-likelihood on `data/ebola/20` in point mode
  5. Validate with `just lint-rs`, `just test-rs`, `just generated-check` (new CLI flag), and `just r treetime timetree` on `data/ebola/20` with each mode
