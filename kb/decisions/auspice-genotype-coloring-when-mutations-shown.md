# Auspice JSON lists the Genotype coloring only when a branch shows a mutation

Auspice JSON lists the `gt` coloring (title `Genotype`, categorical) only when at least one node of the written tree has a non-empty `branch_attrs.mutations`. A tree without displayed mutations lists no `gt` coloring.

**Type**: Behavior change in v1.

**Status**: Approved by the maintainer team on 2026-10-04.

**v1 location**: `fn auspice_tree()` in [packages/app-output/src/auspice.rs](../../packages/app-output/src/auspice.rs) builds every node first and then lists the colorings in the order date, excluded, trait, genotype.

## v0 behavior

v0 lists `Genotype` in the colorings of every Auspice JSON it writes, together with `Date`, `Excluded` and `Branch Support` ([packages/legacy/treetime/treetime/CLI_io.py#L288-L301](../../packages/legacy/treetime/treetime/CLI_io.py#L288-L301)), whether or not the tree carries mutations.

## Rule

- The coloring follows the data in the file, never the command that wrote it
- `branch_attrs.mutations` holds nucleotide substitutions and amino-acid mutations. Auspice shows nucleotide substitutions only, so a nucleotide insertion or deletion alone does not make a mutation visible and does not add `gt`
- augur uses the same rule: `are_mutations_defined()` adds `gt` when a branch has a non-empty `mutations` entry (`augur/export_v2.py` of [augur](https://github.com/nextstrain/augur), revision `0d287496`)

## Why

- An empty color option shows nothing: Auspice colors every branch with the same missing value
- Before this rule each command decided separately: `ancestral` listed `gt` when any mutation existed, `optimize` always, `prune` and `timetree` when they had a root sequence, `clock` and `mugration` never

## Impact

- Auspice files of `optimize`, `prune` and `timetree` runs without a displayed mutation no longer list `gt`
- An `ancestral` tree whose only nucleotide mutations are insertions or deletions no longer lists `gt`
- Trees with at least one displayed mutation keep `gt`, after the date, excluded and trait colorings
