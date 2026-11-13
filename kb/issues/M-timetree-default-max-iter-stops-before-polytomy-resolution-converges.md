# The default `--max-iter 2` stops the time tree before the polytomy resolution converges

With polytomy resolution on (the default), the refinement loop of `treetime timetree` counts a round as converged only when the round resolves no node and the largest node-time change `max_dt` is below 0.01 y. With the default `--max-iter 2`, every measured run stops at the round limit while round 2 still resolves nodes, so the depth of the resolution in the default output, and the dates that depend on it, are set by the round limit and not by convergence. v0 has the same default (`--max-iter` 2, `packages/legacy/treetime/treetime/argument_parser.py#L298`) and the same stop condition (`packages/legacy/treetime/treetime/treetime.py#L374`); the parameterization is approved and does not change.

## Evidence

Measured with `treetime timetree -v` at commit `1e8fbcd8b`, default settings, release build:

| Dataset | Seeds | Rounds of the default run | Nodes resolved in round 2 | With `--max-iter 10` (seeds 1-5) |
| --- | --- | --- | --- | --- |
| `ebola/100` | 1-20 | 2 (limit) | 1-29 over all datasets | converged after 5-10 rounds, 28-38 nodes added; root 0.02-0.06 y earlier, clock rate 7.44e-4 to 7.64e-4 (default runs: 7.608e-4 to 7.667e-4) |
| `mpox/clade-ii/500` | 1-20 | 2 (limit) | | not converged in 10 rounds, 112-120 nodes added |
| `flu/h3n2/500` | 1-20 | 2 (limit) | | converged after 8-10 rounds; see [M-timetree-root-date-moves-between-refinement-rounds-on-flu-h3n2-500.md](M-timetree-root-date-moves-between-refinement-rounds-on-flu-h3n2-500.md) |

- After 2 rounds, 7-57 polytomies remain, and the count falls in each later round
- Code: the round loop in `packages/treetime/src/timetree/refinement_loop.rs`, the convergence record in `packages/treetime/src/timetree/convergence/optimizer.rs`

> [!IMPORTANT]
> **Investigation required.** v0 was not run on these datasets. Run v0 with its defaults on `ebola/100` and `mpox/clade-ii/500` and check whether its tree also stays unconverged after 2 rounds; only a difference from v0 is a v1 defect. Until then, the issue documents the effect of the default for users: report it in the run log or the docs, so that a user who needs a converged resolution raises `--max-iter`.

## Validation

- v0 and v1 give the same round count and the same convergence state on the measured datasets, or a recorded difference explains the divergence
