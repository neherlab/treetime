# Partition Assembly

A `--partition` entry takes a set of sites of one alignment file ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)). Each partition is reconstructed in its own coordinates: local position i of a partition is the i-th of its columns, in ascending order. The outputs, the indel events and the sparse topology decisions belong to the alignment file, so the results of the partitions of one file are mapped back to file coordinates. This page describes that mapping, the assembly of a split alignment.

The owner is `struct AlignmentLayout` ([packages/treetime/src/partition/assembly.rs#L59-L62](../../packages/treetime/src/partition/assembly.rs#L59-L62)), the one owner of the mapping between partition and alignment coordinates. It has two forms (private `enum AlignmentColumns`):

- **Whole** (`label`, `length`): the alignment has one partition over all its columns. The layout holds no per-site data, so a run without `--partition` entries pays nothing per column
- **Split** (`partitions`, `sites`): for every column of the alignment, the partition that takes it and the local position in that partition, and for every partition its columns

`fn alignment_groups()` ([packages/treetime/src/partition/assembly.rs#L11-L40](../../packages/treetime/src/partition/assembly.rs#L11-L40)) takes one pair of an optional partition layout and a partition length per partition, groups the partitions by alignment and builds one `AlignmentLayout` per alignment: a whole layout for a single partition without a layout or with a layout over every column (`fn AlignmentLayout::new()`), otherwise a split layout from the layouts of the partitions. `fn PartitionSet::alignment_groups()` builds the groups of the partitions of a set. Each sparse topology pass builds the groups once and passes them to every edge operation of the pass.

## Whole alignments pass through

An alignment with one partition over all its columns (`fn AlignmentLayout::is_whole()`) has the whole layout and passes the results of that partition through unchanged: sequences, mutations, indel events and ambiguity codes. Assembly then maps no position and builds no list of columns. Runs without partition entries therefore keep their events and their output order.

## Sequences, ambiguity and substitutions

- **Node sequence** (`fn assemble_sequence()`): the sequence of the alignment's length, with each partition's node sequence written to its columns
- **Ambiguity** (`fn assemble_ambiguity()`): the union of the ambiguous sites of the partitions (`type AmbiguousSites`, site position to ambiguity code), each position mapped through the columns of its partition. Two partitions can report an ambiguity at the same local position; each lands on its own alignment site
- **Substitutions** (`fn assemble_substitutions()`): the position of each substitution mapped through the columns of its partition, then sorted by position

Assembly checks that it receives one value per partition and that every position fits its partition; a mismatch is an internal error.

## Indels

An indel is a half-open range of alignment positions with the characters it removes or adds (`struct InDel`, [packages/treetime/src/seq/indel.rs#L315-L320](../../packages/treetime/src/seq/indel.rs#L315-L320)). In a split alignment, each partition stores the pieces of an indel that fall in its columns, in local coordinates. Assembly (`fn assemble_indels()`, [packages/treetime/src/partition/assembly.rs#L178-L210](../../packages/treetime/src/partition/assembly.rs#L178-L210)) maps the pieces of an edge to alignment positions and joins them into events. The pieces come from two sources, and each needs its own join rule.

### Pieces of the marginal passes: maximal runs

The marginal passes classify gaps position by position and join consecutive positions into one event (`fn resolve_indels_backward()` and `fn resolve_indels_forward()` in [packages/treetime/src/seq/indel.rs](../../packages/treetime/src/seq/indel.rs)). Their indels never hold two adjacent events of one kind. Assembly maps the positions of such pieces to alignment positions, splits them into blocks of consecutive alignment positions, and joins the blocks of one kind (deletion or insertion) into maximal runs (`fn join_runs()`, [packages/treetime/src/partition/assembly.rs#L454](../../packages/treetime/src/partition/assembly.rs#L454)). On the same tree, an unsplit alignment and the same alignment split into partitions then give the same indels. For example, a deletion of alignment sites 10 to 12 in an alignment split into codon positions is three one-site pieces, one per partition, and assembly joins them into the one deletion `[10,13)`.

### Pieces split from a composed indel: the identity of the indel

The sparse topology operations compose the indels of two edges with `fn compose_indels()` ([packages/treetime/src/seq/indel.rs#L9](../../packages/treetime/src/seq/indel.rs#L9)). Composition merges adjacent deletions, but keeps adjacent insertions as separate events. A maximal-run join would turn a parent insertion at `[0,1)` and a child insertion at `[1,2)` into one event, and change the event counts, the indel rates and later topology decisions.

`fn split_indels()` ([packages/treetime/src/partition/assembly.rs#L233](../../packages/treetime/src/partition/assembly.rs#L233)) splits a list of alignment indels into the pieces of each partition: each partition gets the positions of the indel in its columns, joined into runs of consecutive local positions. Every piece records the identity of its alignment indel in `InDel.split_from` (`struct SplitFrom { index, range }`, [packages/treetime/src/seq/indel.rs#L354-L357](../../packages/treetime/src/seq/indel.rs#L354-L357)): the index of the indel in the split list and its alignment range. Only `split_indels()` sets the field, and only for an alignment with two or more partitions. Assembly joins the pieces of one identity into exactly that indel (`fn join_split_pieces()`), never joins pieces of two identities, and checks that the pieces cover the recorded range.

An edge holds pieces of one source only. Pieces of the marginal passes next to split pieces on one edge are an internal error, because every operation that changes the indels of an edge replaces all its pieces with `split_indels()` of the new alignment list.

### Invariant

For every range-sorted list `L` of alignment indels of an edge, in every layout:

```text
assemble(split(L)) == L
```

The unit and property tests of `packages/treetime/src/partition/__tests__/test_assembly.rs` check it in interleaved and block layouts, including a deletion and an insertion across partition borders, and the two adjacent insertions `[0,1)` and `[1,2)`, which stay two events.

## Mutation order of an edge

The writers (Newick comments, Auspice JSON, UShER MAT) emit the mutations of an edge in stored order, so assembly keeps the order of the unsplit alignment (`fn assemble_edge_mutations()`, combined with `fn combine_edge_mutations()` of [packages/treetime/src/seq/mutation.rs#L174](../../packages/treetime/src/seq/mutation.rs#L174)):

1. The substitutions, by position
2. The indels, in the order in which the unsplit alignment stores them on the same tree. Indels of the marginal passes come as all deletions, then all insertions, each by range, the order of `fn resolve_indels_forward()`. Indels that a collapse, merge or reroot composed come in the range order that `fn compose_indels()` returns, with the kinds mixed; assembly restores it from `SplitFrom.index`

An edge with an insertion before a deletion in alignment coordinates therefore lists the deletion first when its indels come from the marginal passes.

## Streaming

Ancestral and timetree emit node sequences per alignment (`fn stream_alignment_sequences()`, [packages/treetime/src/partition/alignment_sequences.rs#L18](../../packages/treetime/src/partition/alignment_sequences.rs#L18), the one owner of streaming per alignment). A whole alignment keeps the single pass in which the partition streams its sequences and computes its mutations. A split alignment computes the mutations of each partition without a sink, assembles them, and emits the assembled node sequences and ambiguity codes in preorder, so the run never holds the sequences of all nodes at once.

## Indel counts

The per-edge indel count of the branch-length objective is the number of indel events of each alignment, summed over the alignments: the partition's own events for a whole alignment, the assembled events for a split alignment (`fn gather_alignment_edge_indel_counts()`, [packages/treetime/src/optimize/gather.rs#L98-L120](../../packages/treetime/src/optimize/gather.rs#L98-L120)). The indel rate of an alignment divides its event count by the total branch length ([kb/algo/indel-models.md](indel-models.md)).

## Sparse topology operations

The sparse topology operations work on the assembled indels of each alignment, through the edge-indel owner of sparse reconstructions [packages/treetime/src/partition/marginal/sparse/alignment_indels.rs](../../packages/treetime/src/partition/marginal/sparse/alignment_indels.rs) (assemble the pieces of an edge, compose two edges, store a list of alignment indels as split pieces):

- **Collapse, edge merge and reroot**: assemble the indels of both edges, compose them in alignment coordinates, and store the result as split pieces on the new edge
- **Reversion hoist**: tests overlap and adjacency on the assembled indels of the parent and child edges, so the tests use alignment positions
- **Sibling merge** (`fn build_mutation_index()` in [packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L93](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L93)): keys indels by alignment number and assembled indel, so an indel across several partitions counts once in the shared count; substitutions keep their per-partition keys. A merge moves shared indels by splitting the assembled lists of the edges again, never by matching stored pieces. The child branch lengths of a merged group count each remaining assembled indel of an alignment once

## Related

- [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)
- [kb/algo/indel-models.md](indel-models.md)
