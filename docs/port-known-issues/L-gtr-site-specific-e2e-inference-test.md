# Site-specific GTR inference lacks end-to-end test from real tree data

The site-specific GTR inference is validated against v0 using synthetic mutation counts generated from the model's own equations. An end-to-end test that derives per-site mutation counts (`n_ija`, `T_ia`, `root_state`) from a real tree and alignment via the partition system is missing.

## Current test coverage

- Self-consistency test: synthetic counts from model equation, solver recovery at 1e-3
- Golden-master test: same synthetic counts compared between v0 and v1 solvers at 1e-3
- Property tests: column-stochastic, semigroup, equilibrium, etc. at 1e-8 to 1e-10

## What the end-to-end test would catch

- Parent/child orientation mistakes in the mutation count extraction path
- Interaction between site-specific GTR and the branch joint distribution computation
- Midpoint counting formula compatibility between dense partition edge messages and site-specific inference input

## Dependency

Requires partition integration (`L-gtr-site-specific-partition-integration.md`). The mutation count extraction functions (`get_branch_mutation_matrix`, `accumulate_mutation_counts`) currently produce site-summed `MutationCounts`, not per-site `MutationCountsSiteSpecific`.
