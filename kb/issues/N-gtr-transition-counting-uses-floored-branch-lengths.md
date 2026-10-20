# GTR transition counting uses floored branch lengths

The dense and discrete transition counts apply the branch-length floor of the marginal passes before they evaluate `expQt` and accumulate dwell times: `fn count_transitions_dense()` calls `effective_branch_length()` [packages/treetime/src/partition/marginal/shared/data.rs#L31](../../packages/treetime/src/partition/marginal/shared/data.rs#L31) [packages/treetime/src/partition/marginal/shared/data.rs#L64-L68](../../packages/treetime/src/partition/marginal/shared/data.rs#L64-L68). The sparse count clamps the same way [packages/treetime/src/partition/marginal/sparse/count.rs#L26-L33](../../packages/treetime/src/partition/marginal/sparse/count.rs#L26-L33).

GTR estimation therefore sees floored branch lengths instead of raw values. Zero-length branches contribute dwell time `T_i` at the floor length, which can bias the rate matrix toward short-branch statistics. No bias is measured.

> [!IMPORTANT]
> **Decision required.** The behavior reproduces v0, which floors the branch length for GTR mutation counting through `_branch_length_to_gtr()` [packages/legacy/treetime/treetime/treeanc.py#L752-L760](../../packages/legacy/treetime/treetime/treeanc.py#L752-L760), used at [packages/legacy/treetime/treetime/treeanc.py#L1106](../../packages/legacy/treetime/treetime/treeanc.py#L1106). The options are:
>
> - Keep v0 parity and close this issue
> - Count transitions with raw branch lengths, handle zero-length branches explicitly, and record the divergence from v0 in `kb/decisions/`
