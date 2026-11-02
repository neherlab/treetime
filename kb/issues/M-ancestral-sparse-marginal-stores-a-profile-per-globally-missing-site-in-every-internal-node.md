# Sparse marginal reconstruction stores a profile per globally missing site in every internal node

A column without data in any tip (globally missing: `N`, `-`, or another non-character in every sequence) has the equilibrium frequencies $\pi$ as its posterior at every node, because $\pi P(t) = \pi$ for every branch length $t$. The sparse marginal passes still store one per-site profile for each such column in every internal node, so memory and run time grow with $N \cdot K$ for $N$ internal nodes and $K$ globally missing columns.

## Mechanism

- **Root**: the backward pass gives the root the intersection of the no-data ranges of its children (`fn node()`, [packages/treetime/src/partition/marginal/sparse/backward.rs#L85](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L85)); a leaf contributes its gaps and unknown ranges ([#L81](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L81)). `fn with_equilibrium_at()` ([packages/treetime/src/partition/marginal/sparse/backward.rs#L231-L245](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L231-L245)), called at [#L108](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L108), expands every root no-data range into single `VarPosMap` entries holding $\pi$, one per position, ranges of gaps included
- **Forward pass**: `fn compute_msg_to_child()` collects an outside profile for every position of the child's no-data ranges ([packages/treetime/src/partition/marginal/sparse/forward.rs#L272-L291](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L272-L291)). For an internal child the flag `outside_at_missing` is always set ([#L136](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L136)). Where the parent's sequence holds a gap or another non-character, the `else if` branch at [#L287](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L287) copies the parent's per-site entry. `fn edge()` then propagates these entries along the branch and merges them into the profile of the internal child ([#L176-L190](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L176-L190)). A globally missing column lies in the no-data ranges of every node, so every internal node receives and propagates one entry per column

Commit `430a86373` ("infer internal states without descendant data from the whole tree") added the root expansion and the outside entries at the no-data sites of internal nodes, so that an internal node without descendant data takes the most likely state of its posterior, as v0 does ([kb/decisions/ancestral-marginal-tip-reconstruction-and-imputation.md](../decisions/ancestral-marginal-tip-reconstruction-and-imputation.md)). For a clade without data below a node, the posterior depends on the path from the rest of the tree, and a per-site entry is needed. For a globally missing column, the posterior is $\pi$ at every node, and its most likely state is the same everywhere: on `sc2/2844` with masked columns, all 1951 internal nodes write `T`, the most likely state of $\pi$, where `d45ee9861` wrote `N`. Where the tips hold gaps, the internal nodes inherit the gap, so the propagated profile is not used for the output at all: at column 95 of `lassa/L/500` (gaps in 413 tips, `N` in one), all 412 internal nodes write `-`.

Commit `2a41a522e` removed the shared profiles of sites without data below a child, which grouped such sites by their pattern of children without data, and recorded speed and memory on alignments with long runs of `N` as the cost.

## Measurement

Release build, `treetime ancestral --method-anc=marginal --representation=sparse -j 4` on omega, peak resident memory from `/usr/bin/time`. "Before" is `430a86373^` (`d45ee9861`), "after" is `4535dcd01`. The masked inputs replace columns of `data/sc2/2844/aln.fasta.xz` with `N` in every sequence. Of the example datasets in `data/`, only `ebola`, `lassa/L`, and `mpox/clade-ii/20` have globally missing columns.

| Input | Internal nodes | Globally missing columns | Before: time, memory | After: time, memory |
| --- | ---: | ---: | --- | --- |
| `sc2/2844` | 1951 | 0 | 1.51 s, 0.61 GB | 1.54 s, 0.58 GB |
| `sc2/2844`, columns 10 001-20 000 masked | 1951 | 10 000 | 1.86 s, 0.58 GB | 4.90 s, 1.80 GB |
| `sc2/2844`, columns 5 001-25 000 masked | 1951 | 20 000 | 1.81 s, 0.48 GB | 7.53 s, 2.99 GB |
| `lassa/L/500` | 412 | 2 334 | 0.80 s, 0.22 GB | 1.09 s, 0.28 GB |

The extra memory grows linearly with the number of masked columns: about 64 bytes per internal node and globally missing column (1.22 GB for $1951 \times 10\,000$ entries, 2.51 GB for $1951 \times 20\,000$).

## Fix direction

Represent the globally missing columns of a partition once, as ranges with the shared profile $\pi$, instead of per-site entries: the root keeps the ranges, the forward pass skips them, and the reconstruction writes the most likely state of $\pi$ (or the gap that the internal nodes inherit where the tips hold gaps) at every internal node. Clade-level no-data sites keep their per-site entries. Compare the states at these columns with the dense backend and with the outputs of `430a86373` on the masked `sc2/2844` input.

## Related issues

- [M-ancestral-sparse-stores-a-full-sequence-per-node.md](M-ancestral-sparse-stores-a-full-sequence-per-node.md): the per-node sequences of the sparse backend, an $O(N \cdot L)$ cost of its own that this issue does not include
- [N-ancestral-sparse-profiles-allocate-for-large-alphabets.md](N-ancestral-sparse-profiles-allocate-for-large-alphabets.md): each per-site entry allocates on the heap with an amino-acid alphabet, which multiplies the cost of this issue for protein partitions
