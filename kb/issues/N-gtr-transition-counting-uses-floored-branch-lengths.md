# GTR transition counting uses floored branch lengths

The dense and discrete transition counts apply the branch-length floor of the marginal passes before they evaluate `expQt` and accumulate dwell times: `fn count_transitions_dense()` calls `effective_branch_length()` [packages/treetime/src/partition/marginal/shared/data.rs#L29](../../packages/treetime/src/partition/marginal/shared/data.rs#L29) [packages/treetime/src/partition/marginal/shared/data.rs#L63-L65](../../packages/treetime/src/partition/marginal/shared/data.rs#L63-L65).

GTR estimation therefore sees floored branch lengths instead of raw values. Zero-length branches contribute dwell time `T_i` at the floor length, which can bias the rate matrix toward short-branch statistics. No bias is measured.

> [!IMPORTANT]
> **Decision required.** The behavior reproduces v0, which floors the branch length for GTR mutation counting through `_branch_length_to_gtr()` [packages/legacy/treetime/treetime/treeanc.py#L752-L760](../../packages/legacy/treetime/treetime/treeanc.py#L752-L760), used at [packages/legacy/treetime/treetime/treeanc.py#L1106](../../packages/legacy/treetime/treetime/treeanc.py#L1106). The options are:
>
> - Keep v0 parity and close this issue
> - Count transitions with raw branch lengths, handle zero-length branches explicitly, and record the divergence from v0 in `kb/decisions/`

## Proposed implementation plan

`fn count_transitions_dense()` floors at [packages/treetime/src/partition/marginal/shared/data.rs#L29](../../packages/treetime/src/partition/marginal/shared/data.rs#L29) through `fn DenseInputs.effective_branch_length()` (`raw.max(self.min_branch_length)`), and v0 floors the same way at the cited lines. The plan stops at the open decision; one short plan per option:

- **Keep v0 parity**: delete this issue; no code change.
- **Raw branch lengths**: in `fn count_transitions_dense()`, use the raw length for `Ti` and for the transition matrix, and give zero-length branches zero dwell time and no transitions (their child state equals the parent state in the expectation); write the divergence into a new `kb/decisions/` file with approval. Tests: unit test on a 3-tip tree with one zero-length branch where `Ti` equals the hand-computed sum of raw lengths (oracle: manual calculation), and a simulation measurement of the `W` and `pi` error with and without the floor (oracle: the true simulated model) to show the bias that this issue assumes.

Validation: `just lint-rs`, `just test-rs`.
