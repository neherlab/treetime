# Restrict amino-acid mutation output to substitutions

## Motivation

TreeTime v1 emits amino-acid mutations as substitutions plus range-based insertion and deletion events. The nucleotide track was recently made substitution-only so that the Auspice `mutations.nuc` list matches the augur node-data contract, but the amino-acid track still carries indel tokens. This proposal records the question of whether the amino-acid track should likewise be restricted to substitutions, and the tension that decision has with reference parity.

The two tracks are produced by different code:

- Nucleotide node-data `muts` are substitution-only: `fn build_augur_node_data_json()` maps each `Sub` to its string and never emits indels [`packages/treetime/src/commands/ancestral/augur_node_data.rs#L71-L77`](../../packages/treetime/src/commands/ancestral/augur_node_data.rs#L71-L77).
- Amino-acid node-data chains substitutions with range-based indels and serializes both: [`packages/treetime/src/commands/ancestral/aa_node_data.rs#L270-L276`](../../packages/treetime/src/commands/ancestral/aa_node_data.rs#L270-L276). The substitution diff excludes any position where either side is a gap, so gaps reach the output only through the indel track: `fn is_reportable_sub()` [`packages/treetime/src/commands/ancestral/aa_node_data.rs#L326-L329`](../../packages/treetime/src/commands/ancestral/aa_node_data.rs#L326-L329).
- The Auspice writer groups whatever mutations it is given, per track: `fn group_mutations()` [`packages/treetime/src/commands/shared/tree_output.rs#L1473-L1502`](../../packages/treetime/src/commands/shared/tree_output.rs#L1473-L1502). It now drops indels for the nucleotide track only; amino-acid tracks pass through `fn mutation_event_strings()` [`packages/treetime/src/seq/mutation.rs#L57`](../../packages/treetime/src/seq/mutation.rs#L57), which expands a range-based indel into one positional token per residue (`{state}{pos}-` for a deletion, `-{pos}{state}` for an insertion).

## Reference behavior (augur 34.0.0, 2026-07-07)

augur does not model insertions or deletions as range events on either track. It compares aligned sequences position by position and writes each difference as a single `<from><pos><to>` string, allowing `-` on either side as an ordinary aligned character.

Nucleotide muts come from the TreeTime substitution list, one string per differing site [[src](https://github.com/nextstrain/augur/blob/d8e38736037ba9474a809f9a5a63bc2b279d2407/augur/ancestral.py#L218)] inside `def collect_mutations()` [[src](https://github.com/nextstrain/augur/blob/d8e38736037ba9474a809f9a5a63bc2b279d2407/augur/ancestral.py#L185)].

Amino-acid muts come from a direct column comparison of the translated reference against each node's translation, with no gap filtering [[src](https://github.com/nextstrain/augur/blob/d8e38736037ba9474a809f9a5a63bc2b279d2407/augur/translate.py#L305)]:

```python
muts = [construct_mut(a, int(pos+1), d) for pos, (a,d) in enumerate(zip(ref, aln[n.name])) if a!=d]
```

`def construct_mut()` simply concatenates the three fields [[src](https://github.com/nextstrain/augur/blob/d8e38736037ba9474a809f9a5a63bc2b279d2407/augur/translate.py#L206)], and translation renders in-frame gaps as `-` [[src](https://github.com/nextstrain/augur/blob/d8e38736037ba9474a809f9a5a63bc2b279d2407/augur/translate.py#L67)]. A deleted residue therefore appears as `K100-`, and a gap-to-residue change as `-100A`, at fixed alignment-column coordinates.

augur's Auspice export copies the node-data muts verbatim into the tree; `def create_branch_mutations()` [[src](https://github.com/nextstrain/augur/blob/d8e38736037ba9474a809f9a5a63bc2b279d2407/augur/export_v2.py#L1062)] assigns `branch_attrs[...]["mutations"]["nuc"] = node_info["muts"]` [[src](https://github.com/nextstrain/augur/blob/d8e38736037ba9474a809f9a5a63bc2b279d2407/augur/export_v2.py#L1070)]. Node-data and Auspice mutation lists are the same object by construction and cannot disagree.

## The divergence question

augur amino-acid muts are **not** substitution-only: a deletion is present as a `K100-` gap-substitution. TreeTime's token format matches (a range deletion expands to the same `K100-` strings), so for deletions the two representations are string-compatible when coordinates align. The mismatch risk is confined to:

1. **Insertions.** augur has no insertion concept on a fixed-length alignment; an inserted residue relative to the reference is not a separate column. TreeTime's range-based insertions serialize as `-{pos}{state}`, whose position is a reference coordinate that may not correspond to any augur alignment column. Amino-acid sequences are translated from a fixed nucleotide MSA, so true insertions are not expected to arise, but the code path can emit them.
2. **Coordinate correspondence.** Range indel coordinates are 1-based over the reference; augur coordinates are 1-based over the translated alignment column. Whether these coincide for every CDS has not been verified against augur output.

Restricting the amino-acid track to substitutions removes both risks, but drops the deletion tokens that augur _does_ emit as gap-substitutions. It trades one divergence (possibly wrong indel coordinates) for another (missing deletions).

## Design axes

Two independent decisions.

### A1. Amino-acid indel inclusion

- **A1a. Substitutions only.** Drop amino-acid indels entirely, mirroring the nucleotide track. Simplest and internally consistent, avoids any coordinate mismatch, but omits deletions that augur represents.
- **A1b. Keep range-based indels (current).** Retains deletion and insertion information but must prove coordinate correspondence with augur before it can claim parity.
- **A1c. Emit amino-acid gaps as positional gap-substitutions.** Reproduce augur's model exactly: derive `K100-` / `-100A` from a per-column comparison of translated sequences rather than from range events. Closest reference parity; requires routing amino-acid gaps through the substitution path instead of the indel track.

### A2. Scope of the change

- **A2a. Both node-data and Auspice.** Keep the two outputs identical (augur's invariant). Any restriction applies at the amino-acid node-data source so every consumer agrees.
- **A2b. Auspice only.** Filter in `fn group_mutations()`, leaving amino-acid node-data unchanged. This reintroduces exactly the node-data-vs-Auspice disagreement that the nucleotide fix removed and is not recommended.

## Recommendation (open)

Decision required from the team. The narrow request that motivated this proposal is A1a (amino-acid substitutions only), applied per A2a so node-data and Auspice stay identical. That is defensible as a consistency choice and is safe against the insertion-coordinate risk, but it is a deliberate divergence from augur, which surfaces amino-acid deletions. A1c is the only option that matches augur exactly. This proposal does not select an option; exact parity with augur remains the default until the team approves a divergence.

## Impact

- CLI outputs affected: `ancestral`, `optimize`, `timetree` amino-acid node-data and Auspice v2 `mutations` for CDS tracks. Datasets with amino-acid indels (translated from indel-bearing nucleotide alignments) change; substitution-only datasets are unaffected.
- No nucleotide-track change; the nucleotide substitution-only contract already holds.

## Validation plan

- Golden-master comparison of amino-acid `muts` against `augur translate` node-data for an indel-bearing CDS, per option, to measure the actual divergence.
- Unit coverage in `fn group_mutations()` and the amino-acid node-data builder asserting the chosen per-track policy.
- Confirm node-data and Auspice amino-acid lists remain identical under A2a.

## Non-goals

- No change to the nucleotide track.
- No change to PhyloXML or UShER MAT indel encoding, which have their own format contracts.

## Related

- [kb/issues/N-amino-acid-mutation-indel-representation-undecided.md](../issues/N-amino-acid-mutation-indel-representation-undecided.md)
- [phyloxml-treetime-property-namespace.md](phyloxml-treetime-property-namespace.md)
