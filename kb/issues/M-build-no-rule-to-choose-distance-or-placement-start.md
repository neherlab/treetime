# Tree builder has no rule to choose between a distance start and placement

The prototype tree builder can start from sequential placement or from neighbor joining (`InitialTree` in [`packages/treetime/src/build/pipeline.rs#L64`](../../packages/treetime/src/build/pipeline.rs#L64)). The user picks the start. No rule chooses it from the data.

## Evidence

Fitch score relative to the reference tree, after rebuilds of radius 3 ([kb/reports/feat-tree-infer-prototype.md](../reports/feat-tree-infer-prototype.md)):

| Dataset | Mean p-distance to consensus | Placement, random order | Neighbor joining |
| --- | --- | --- | --- |
| `rsv/a/2000` | 0.014 | −0.01% | +0.14% |
| `flu/h3n2/500` | 0.022 | −0.07% | −0.07% |
| `dengue/2000` | 0.049 | +0.05% | +0.12% |
| `lassa/L/500` | 0.121 | +3.09% | +0.96% |

The crossover lies between 0.049 and 0.121, with only one dataset above it.

## Impact

Without a rule, a default command either leaves Lassa-like data about 3% above the reference tree or pays the cubic cost of neighbor joining on low-divergence data.

## Fix

Simulate datasets with $\mu$ = 0.02 to 0.03 and gamma rates (`examples/simulate_dataset.rs`), add further divergent real datasets, and place a threshold on the mean p-distance. The options are listed in [kb/proposals/native-parsimony-tree-builder.md](../proposals/native-parsimony-tree-builder.md), axis B2. The choice changes the inferred topology, so it needs a decision file.
