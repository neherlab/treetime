# Auspice JSON without another default colors the tree by its most recurrent site

When an Auspice JSON has no default coloring from its filters (no `bad_branch` from dates and no trait), `display_defaults.color_by` names a genotype coloring `gt-nuc_<position>`. The position is the site with substitutions between determined states (A, C, G, T for nucleotides) on the most branches. Ties go to the lowest position. A tree without such a substitution gets no default coloring.

**Type**: Behavior change in v1.

**Status**: Approved by the maintainer on 2026-10-06.

**v1 location**: `fn auspice_data()` and `fn genotype_color_by()` in [packages/app-output/src/auspice.rs](../../packages/app-output/src/auspice.rs). The ranking is `fn sites_by_branch_count()` in [packages/treetime/src/homoplasy/site_branches.rs](../../packages/treetime/src/homoplasy/site_branches.rs), which also ranks the "Sites mutated on several branches" table of the app's ancestral results, so the default site is the first row of that table.

## v0 behavior

v0 writes the `Excluded` coloring (`bad_branch`) into every Auspice JSON and sets `display_defaults.color_by` to `bad_branch` ([packages/legacy/treetime/treetime/CLI_io.py#L288-L303](../../packages/legacy/treetime/treetime/CLI_io.py#L288-L303)). An `ancestral` run marks no branch as excluded, so every branch has the value `No` and the tree has one color.

## Why

- Auspice starts at its built-in default coloring `country`. When the file lacks it, Auspice falls back to `display_defaults.color_by`, then to the first coloring other than `gt`. A genotype coloring needs a position, so Auspice never falls back to `gt`
- v1 `ancestral` and `optimize` files list only the `gt` coloring ([auspice-genotype-coloring-when-mutations-shown.md](auspice-genotype-coloring-when-mutations-shown.md)). Without a default, Auspice logs an error, sets the coloring to `none` and draws a grey tree with an `unknown` legend
- Auspice accepts a genotype value in `display_defaults.color_by` and checks its position against `meta.genome_annotations` (`checkAndCorrectErrorsInState()` in `src/actions/recomputeReduxState.js` of [auspice](https://github.com/nextstrain/auspice), tag `v3.0.0`)
- The site with the most independent substitutions is the most informative single site to show first: it is the strongest homoplasy candidate of the tree

## Impact

- `ancestral`, `optimize` and `homoplasy` files with a substitution between determined states open colored by that site, in the app and in Auspice
- Files with dates or a trait keep the filter default (`bad_branch`, or the trait)
- `prune` files without mutations and files whose only nucleotide changes involve ambiguous characters keep no default and still open grey
- If `ancestral` output adopts the v0 `bad_branch` coloring ([kb/issues/N-ancestral-auspice-json-not-produced.md](../issues/N-ancestral-auspice-json-not-produced.md), axis A2), the filter default takes precedence and this rule no longer applies to `ancestral`
