# Reroot split-optimizer default diverges from v0

Clock and timetree use an $11$-point grid as the default continuous split-position search, while v0 brackets on a grid and refines the selected edge with bounded scalar minimization [packages/legacy/treetime/treetime/treeregression.py](../../packages/legacy/treetime/treetime/treeregression.py). Optimize already selects Brent refinement explicitly. Changing the clock/timetree default would alter numerical output and requires approval independently of code-structure migration.

## Options

- O1. Match v0 by making bounded Brent refinement the clock/timetree default, with the same bracket and stopping contract.
- O2. Retain the $11$-point grid and document an intentional numerical divergence. This requires evidence that the coarser root fraction is scientifically acceptable.

## Recommendation

Use O1 for reference parity after differential tests establish the exact v0 bracket and tolerance. No implementation ticket is ready until this numerical-behavior decision is approved.

## Related issues

- [N-reroot-clock-search-duplicates-generic-module.md](N-reroot-clock-search-duplicates-generic-module.md)
