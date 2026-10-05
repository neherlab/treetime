# Marginal branch-length clamping has unresolved numerical policy and ownership

> [!WARNING]
> **Needs review.** Two claims about the clamp are not verified against the current code. (1) The non-finite failure at zero branch length was not reproduced. The dense backward pass takes the logarithm of each child message [packages/treetime/src/partition/marginal/shared/pass.rs#L106-L108](../../packages/treetime/src/partition/marginal/shared/pass.rs#L106-L108), and at $t=0$ a leaf message keeps the exact zeros of its one-hot profile, so a likely failure is $\ln 0 = -\infty$ in every state when children disagree. That mechanism is not confirmed by a test run. (2) Zero-branch-length tests are said to have been weakened from exact properties (`Ti == 0` with epsilon `1e-15`, off-diagonal mass `< 0.1`) to `is_finite()` checks. The cited file `packages/treetime/src/ancestral/__tests__/test_marginal_dense.rs` no longer exists, and no current test with those assertions was found. The closest current test checks only finiteness, sign, and a `1e-2` neighbourhood of $\ln \pi_A$ at $t = 10^{-8}$ [packages/treetime/src/ancestral/__tests__/test_marginal_branch_length/test_marginal_branch_length_equilibrium.rs#L88-L106](../../packages/treetime/src/ancestral/__tests__/test_marginal_branch_length/test_marginal_branch_length_equilibrium.rs#L88-L106).

> [!IMPORTANT]
> **Decision required.** The floor matches v0, which raises every marginal branch length to `MIN_BRANCH_LENGTH * one_mutation` with `MIN_BRANCH_LENGTH = 1e-3` [packages/legacy/treetime/treetime/treeanc.py#L752-L760](../../packages/legacy/treetime/treetime/treeanc.py#L752-L760). Keeping the clamp preserves v0 parity. Exact-zero handling ($e^{Q \cdot 0} = I$ through a log-safe path) or rejection of zero-length branches would diverge from v0 and needs approval with parity and numerical validation.

The marginal passes raise every branch below `MIN_BRANCH_LENGTH_FRACTION / sequence_length` to that threshold before message propagation. The same floor is implemented in three places:

- `fn fix_branch_length()` [packages/treetime/src/hacks/fix_branch_length.rs#L5-L8](../../packages/treetime/src/hacks/fix_branch_length.rs#L5-L8), called only by the sparse passes [packages/treetime/src/partition/marginal/sparse/forward.rs#L108](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L108) [packages/treetime/src/partition/marginal/sparse/backward.rs#L139](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L139)
- the dense and discrete passes, through `DenseInputs::min_branch_length` [packages/treetime/src/partition/marginal/dense/partition.rs#L64](../../packages/treetime/src/partition/marginal/dense/partition.rs#L64), applied inline in [packages/treetime/src/partition/marginal/shared/pass.rs#L132](../../packages/treetime/src/partition/marginal/shared/pass.rs#L132) and [packages/treetime/src/partition/marginal/shared/pass.rs#L252](../../packages/treetime/src/partition/marginal/shared/pass.rs#L252)
- sparse transition counting [packages/treetime/src/partition/marginal/sparse/count.rs#L29-L36](../../packages/treetime/src/partition/marginal/sparse/count.rs#L29-L36)

The name `fix_branch_length` does not reveal that it clamps, and the `hacks` namespace does not identify the owning algorithm. Mathematically, $e^{Q \cdot 0}=I$ is valid, so changing the clamp is a numerical and scientific decision rather than architecture cleanup.

## Separable work

- Preserve behavior while giving the floor one owner beside marginal inference: a private operation whose name states the floor policy, used by the dense, discrete, and sparse passes. Delete the `hacks` module afterwards
- Decide clamping, exact-zero handling, or rejection only with parity and numerical validation (see the decision block above)
