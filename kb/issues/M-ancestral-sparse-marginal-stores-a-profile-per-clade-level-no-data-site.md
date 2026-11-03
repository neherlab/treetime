# Sparse marginal reconstruction stores a profile per clade-level no-data site in every internal node of the clade

A clade-level no-data site of an internal node is a position where the subtree of the node has no data while other tips of the tree have data there. Its posterior at the node is the message from the rest of the tree, which depends on the path to the node, so the sparse forward pass stores one per-site entry for each such site in every internal node of the clade. Memory and run time grow with the number of internal nodes times the length of their no-data ranges, for example with long runs of `N` in a clade of tips sequenced with the same primer scheme.

Columns without data in any tip are not part of this cost: their posterior is $\pi$ at every node, and the sparse passes store their ranges once per update with $\pi$ (`struct GloballyMissing` in [packages/treetime/src/partition/marginal/sparse/globally_missing.rs](../../packages/treetime/src/partition/marginal/sparse/globally_missing.rs)), as [kb/algo/ancestral.md](../algo/ancestral.md) describes under "Globally missing columns".

## Mechanism

`fn compute_msg_to_child()` ([packages/treetime/src/partition/marginal/sparse/forward.rs#L284-L305](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L284-L305)) collects the undivided parent posterior at every position of the child's no-data range outside the globally missing columns. For an internal child the flag `outside_at_missing` is always set ([#L138](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L138)), and `fn edge()` propagates these entries along the branch and merges them into the profile of the child ([#L179-L194](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L179-L194)). Commit `430a86373` added these entries so that an internal node without descendant data takes the most likely state of its posterior, as v0 does ([kb/decisions/ancestral-marginal-tip-reconstruction-and-imputation.md](../decisions/ancestral-marginal-tip-reconstruction-and-imputation.md)).

Commit `2a41a522e` removed the shared profiles of sites without data below a child, which grouped such sites by their pattern of children without data, and recorded speed and memory on alignments with long runs of `N` as the cost.

## Measurement

Release builds, `treetime ancestral --method-anc=marginal --representation=sparse -j 4` on omega, peak resident memory from `/usr/bin/time`. "Before" is `430a86373^` (`d45ee9861`), the commit without outside entries at no-data sites. A temporary counter at `62465d626` counted the per-site entries in the profiles of internal nodes after the forward pass and split off the entries at globally missing columns, which the current representation no longer stores; the clade-level and the other entries are the same in both representations. "Memory now" is the build with the shared representation of the globally missing columns. The masked input replaces columns 10 001 to 20 000 of `data/sc2/2844/aln.fasta.xz` with `N` in every sequence.

| Input | Internal nodes | Clade-level entries | Other entries | Clade-level share | Memory before | Memory now | Ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `lassa/L/500` | 412 | 392 608 | 290 909 | 57% | 214 MB | 233 MB | 1.09 |
| `ebola/362` | 359 | 102 333 | 1 605 | 98% | 80 MB | 85 MB | 1.05 |
| `mpox/clade-ii/2000` | 1987 | 3 377 769 | 60 686 | 98% | 3429 MB | 3297 MB | 0.96 |
| `sc2/5100` | 5099 | 1 314 464 | 32 430 | 98% | 1606 MB | 1525 MB | 0.95 |
| `sc2/2844`, masked | 1951 | 318 522 | 6 776 | 98% | 579 MB | 559 MB | 0.97 |

Only `lassa/L/500` and `ebola/362` use more memory than before, by 9% and 5%. On the other inputs, later changes of the sparse backend more than offset the clade-level entries.

> [!IMPORTANT]
> **Decision required.** The clade-level entries cost up to 9% of peak memory on the example datasets. The options:
>
> - Keep one per-site entry per clade-level no-data site, and close this issue as an accepted cost
> - Share one profile among the sites of a clade-level no-data run, grouped by the pattern of children without data. This reverses `2a41a522e`, which removed such shared profiles to keep one kind of shared profile in the sparse passes
>
> Evidence: the table above. A run of `N` in a large clade multiplies the cost by the number of internal nodes of the clade; no example dataset shows a case above 9%.

## Related issues

- [M-ancestral-sparse-stores-a-full-sequence-per-node.md](M-ancestral-sparse-stores-a-full-sequence-per-node.md): the per-node sequences of the sparse backend, an $O(N \cdot L)$ cost of its own that this issue does not include
- [N-ancestral-sparse-profiles-allocate-for-large-alphabets.md](N-ancestral-sparse-profiles-allocate-for-large-alphabets.md): each per-site entry allocates on the heap with an amino-acid alphabet, which multiplies the cost of this issue for protein partitions
