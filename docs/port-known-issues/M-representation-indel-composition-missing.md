# Indel composition missing during edge merge

When two consecutive edges are merged (edge collapse, reroot edge merge), substitutions are properly composed via `fn compose_substitutions()` ([packages/treetime/src/seq/mutation.rs#L72](../../packages/treetime/src/seq/mutation.rs#L72)) but indels are only concatenated. Overlapping or adjacent indels on consecutive edges produce redundant or contradictory indel lists instead of a single net indel.

See also: [algorithm background and composition specification](../port-algo-inventory/indel-models.md#indel-composition-on-consecutive-edges).

## Data model

`struct InDel` ([packages/treetime/src/seq/indel.rs#L6](../../packages/treetime/src/seq/indel.rs#L6)) stores:

- `range: (usize, usize)` - half-open position range `[start, end)` in the alignment
- `seq: Seq` - the nucleotide content occupying that range (length = `end - start`)
- `deletion: bool` - `true`: parent has `seq` at this range, child has gaps. `false` (insertion): parent has gaps, child has `seq`

Indels are stored as `Vec<InDel>` on `struct SparseEdgePartition` ([packages/treetime/src/representation/payload/sparse.rs#L94](../../packages/treetime/src/representation/payload/sparse.rs#L94)). They are not sorted by range; order depends on the Fitch backward pass iteration over `variable_indel` (BTreeMap, sorted by range key) followed by `range_difference` results.

## Example

Three nodes A (root), B (internal), C (leaf). Sequence length 6. Node B is collapsed, merging edges A-B and B-C into A-C.

```
node A:    ACGTCG    root
node B:    A---CG    edge A-B: del (1,4) seq="CGT"
node C:    A----G    edge B-C: del (4,5) seq="C"
```

**Current behavior** (concatenation): merged edge A-C carries both indels.

```rust
vec![
  InDel { range: (1,4), seq: "CGT", deletion: true },  // from edge A-B
  InDel { range: (4,5), seq: "C",   deletion: true },  // from edge B-C
]
```

The final sequence is correct by accident (the ranges do not overlap here), but `Composition::add_indel()` processes each entry independently. When ranges do overlap (see edge cases below), the same positions get subtracted multiple times, producing wrong nucleotide frequency counts.

**Expected behavior** (composition): a single net indel on the merged edge.

```rust
vec![
  InDel { range: (1,5), seq: "CGTC", deletion: true },  // net A->C
]
```

This parallels how `fn compose_substitutions()` ([packages/treetime/src/seq/mutation.rs#L72](../../packages/treetime/src/seq/mutation.rs#L72)) chains A->G + G->T into net A->T at the same position, and cancels A->G + G->A into no net change.

## Edge cases

### Case 1: overlapping deletions (same range)

```
node A:    ACGTAC    root
node B:    A--TAC    edge A-B: del (1,3) seq="CG"
node C:    A---AC    edge B-C: del (1,4) seq="CGT"  <-- seq is WRONG
```

The child edge's `seq` field claims "CGT" was deleted, but positions 1..3 are already gaps in B. The true content at B is "--T", not "CGT". Composition should emit one deletion with the parent's original content:

```
expected:  del (1,4) seq="CGT"    from node A (the actual pre-deletion content)
```

Composition double-counting: `fn Composition::add_indel()` ([packages/treetime/src/seq/composition.rs#L68](../../packages/treetime/src/seq/composition.rs#L68)) subtracts C,G for the first indel and C,G,T for the second, removing C and G twice.

### Case 2: adjacent deletions (non-overlapping, abutting)

```
node A:    ACGTCG    root
node B:    A--TCG    edge A-B: del (1,3) seq="CG"
node C:    A----G    edge B-C: del (3,5) seq="TC"
```

No overlap, no double-counting. Concatenation is correct but produces two entries where one suffices:

```
expected:  del (1,5) seq="CGTC"   single merged deletion
```

### Case 3: overlapping insertions (same range)

```
node A:    A---AC    root (gaps at 1..4)
node B:    ATTTAC    edge A-B: ins (1,4) seq="TTT"
node C:    AGGGAC    edge B-C: ins (1,4) seq="GGG"
```

Parent inserts TTT, child overwrites with GGG. Net effect from A is a single insertion of the child's content:

```
expected:  ins (1,4) seq="GGG"    child sequence wins
```

### Case 4: deletion then insertion (same range, cancellation)

```
node A:    ACGTAC    root
node B:    A---AC    edge A-B: del (1,4) seq="CGT"
node C:    ACGTAC    edge B-C: ins (1,4) seq="CGT"
```

Deletion followed by insertion of the same content. Net effect is no change (cancellation):

```
expected:  (empty)                indels cancel, like sub A->G + G->A
```

### Case 5: deletion then insertion (same range, different content)

```
node A:    ACGTAC    root
node B:    A---AC    edge A-B: del (1,4) seq="CGT"
node C:    ATTTAC    edge B-C: ins (1,4) seq="TTT"
```

Deletion of "CGT" followed by insertion of "TTT". The range transitions from one content to another. The `InDel` type has no "replacement" variant, so the conservative representation keeps both:

```
expected:  del (1,4) seq="CGT"  + ins (1,4) seq="TTT"
```

### Case 6: partially overlapping deletions

```
node A:    ACGTAC    root
node B:    A--TAC    edge A-B: del (1,3) seq="CG"
node C:    A----C    edge B-C: del (2,5) seq="TAC"  <-- seq is WRONG
```

The child edge claims "TAC" was deleted at 2..5, but position 2 is already a gap in B. Decompose into sub-ranges:

```
overlap 2..3:  already gapped by parent, child redundant
child-only 3..5:  del (3,5) seq="TA"   (actual content at B)
parent-only 1..2:  del (1,2) seq="C"   (passes through)

expected:  del (1,5) seq="CGTA"   merged from parent original content
```

### Case 7: non-overlapping (no interaction)

```
node A:    ACGTACGT    root
node B:    A--TACGT    edge A-B: del (1,3) seq="CG"
node C:    A--TAC-T    edge B-C: del (6,7) seq="G"
```

Ranges do not touch. Both pass through unchanged:

```
expected:  del (1,3) seq="CG"  + del (6,7) seq="G"    (same as concatenation)
```

Cases 1, 2, and 7 are expected to dominate in phylogenetic data. Cases 4 and 5 (deletion reversed by insertion on the next branch) require an indel event to be immediately undone, which is biologically unusual on adjacent branches. Case 6 (partial overlap) requires two indel events at similar but not identical positions on adjacent branches.

## Affected call sites

Substitution composition (the correct pattern to follow):

- `fn compose_substitutions()` ([packages/treetime/src/seq/mutation.rs#L72](../../packages/treetime/src/seq/mutation.rs#L72)) - merge-sort on sorted substitution lists with chaining and cancellation
- `fn chain_fitch_subs()` ([packages/treetime/src/representation/payload/sparse.rs#L141](../../packages/treetime/src/representation/payload/sparse.rs#L141)) - convenience wrapper on `SparseEdgePartition` that calls `compose_substitutions()`

Indel concatenation sites to replace with composition:

- `fn collapse_edge()` ([packages/treetime/src/representation/algo/topology_cleanup/collapse.rs#L33](../../packages/treetime/src/representation/algo/topology_cleanup/collapse.rs#L33)) - indel concatenation at [line 62](../../packages/treetime/src/representation/algo/topology_cleanup/collapse.rs#L62): `removed_edge_data.indels` prepended before `child_edge.indels` via `append`
- `fn PartitionMarginalSparse::apply_reroot()` ([packages/treetime/src/representation/partition/marginal_sparse.rs#L133](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L133)) - indel concatenation at [line 191](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L191): `parent_edge.indels` extended with `child_edge.indels`

A third site handles indels differently and is not a direct composition case:

- `fn merge_sibling_pair()` ([packages/treetime/src/commands/prune/run.rs#L449](../../packages/treetime/src/commands/prune/run.rs#L449)) - when creating a new internal node between two siblings, indels are preserved on the child edges unchanged at [line 574](../../packages/treetime/src/commands/prune/run.rs#L574) (only substitutions are split into shared vs remaining). This site does not merge consecutive edges, so composition does not apply directly. No change needed here unless the design evolves to share indels between siblings.

Indel application sites (consumers of the merged list, unchanged but relevant for understanding semantics):

- `fn reconstruct_map_seq()` ([packages/treetime/src/representation/partition/marginal_sparse.rs#L90](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L90)) - applies indels to reconstruct node sequences at [line 101](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L101)
- Fitch backward pass ([packages/treetime/src/commands/ancestral/fitch.rs#L584](../../packages/treetime/src/commands/ancestral/fitch.rs#L584)) - applies indels during ancestral sequence reconstruction

## Impact

- **Composition double-counting**: `fn Composition::add_indel()` ([packages/treetime/src/seq/composition.rs#L68](../../packages/treetime/src/seq/composition.rs#L68)) adjusts character counts per indel entry. Redundant overlapping indels subtract/add the same positions multiple times, producing wrong nucleotide frequency counts. These counts feed into GTR parameter estimation
- **Sequence reconstruction**: currently correct by accident (gap-fill is idempotent for deletions, and insertions overwrite). If indel application order changes or insertions overlap with deletions at the same range, results become order-dependent
- **Accumulation**: repeated topology operations (collapse, reroot, collapse again) can stack three or more redundant indels at the same range

## Composition semantics

See the edge cases section above for worked examples of each case. Summary of composition rules:

| Case | Parent edge    | Child edge     | Result                                        |
| ---- | -------------- | -------------- | --------------------------------------------- |
| 1    | del at range R | del at range R | one del at R with parent's seq                |
| 2    | del + del      | (adjacent)     | merge into one del spanning both ranges       |
| 3    | ins at range R | ins at range R | one ins at R with child's seq                 |
| 4    | del at range R | ins at R, same | cancel (no net indel)                         |
| 5    | del at range R | ins at R, diff | keep both (del + ins, no replacement variant) |
| 6    | del at range R | del at range S | decompose overlap, merge sub-ranges           |
| 7    | any at range R | any at range S | no interaction, both pass through             |

Where R and S denote overlapping (cases 1, 3-6) or non-overlapping (case 7) ranges.

## Theoretical background

Substitution models in phylogenetics are <a id="gloss-use-1"></a>continuous-time Markov chains <sup>[1](#gloss-1)</sup> (CTMCs) on a finite state space (the alphabet). The transition probability matrix $P(t) = e^{Qt}$ satisfies the <a id="gloss-use-2"></a>semigroup property <sup>[2](#gloss-2)</sup> $P(t+s) = P(t) \cdot P(s)$. This is what makes `fn compose_substitutions()` correct: composing substitutions on consecutive edges is equivalent to computing the net transition over the summed branch length. The composition is exact because the state space is finite and the process is memoryless.

The semigroup property does not extend to indel processes. The state space for indel models is the set of all possible sequences (infinite and variable-length), making the matrix exponential intractable. The <a id="gloss-use-3"></a>TKF91 <sup>[3](#gloss-3)</sup> model (<a id="cite-1"></a>[Thorne et al. 1991](https://doi.org/10.1007/BF02193625) [[1](#ref-1)]) is the only indel model with an exact closed-form solution, and it restricts indels to single residues. Multi-residue indel models (TKF92 <a id="cite-2"></a>[Thorne et al. 1992](https://doi.org/10.1007/BF00163848) [[2](#ref-2)], the general geometric indel model, the <a id="gloss-use-4"></a>Poisson Indel Process <sup>[4](#gloss-4)</sup> <a id="cite-3"></a>[Bouchard-Cote and Jordan 2013](https://www.pnas.org/doi/full/10.1073/pnas.1220450110) [[3](#ref-3)]) require approximations: pair HMMs, transducer composition, or ODE-based approaches (<a id="cite-4"></a>[Holmes 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7768254/) [[4](#ref-4)]; <a id="cite-5"></a>[Redelings et al. 2024](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11385596/) [[5](#ref-5)]).

The key difference from substitutions: indel composition is not a simple position-wise operation. Insertions change sequence length, so downstream positions shift. Deletions can overlap with prior insertions or other deletions in ways that have no analogue in substitution composition (where each position is independent).

For treetime's purposes, the indel representation is simpler than the general case because:

- Indels are stored as alignment-coordinate ranges (positions do not shift)
- The Fitch parsimony pass already resolved which positions are gapped vs present
- Composition only needs to merge the observed indel annotations, not compute transition probabilities

This means treetime's indel composition is a bookkeeping operation on alignment coordinates, not a probabilistic inference step. The semigroup property is not directly applicable, but the analogy with substitution composition (chaining, cancellation, passthrough) holds at the level of annotation merging.

## Gotchas and implementation risks

- **Indel sort order is not guaranteed.** Unlike substitutions (sorted by position, at most one per position), indels can have overlapping ranges and are not sorted. The composition function must either sort inputs as a precondition or handle unsorted input. Sorting by `range.0` with ties broken by `range.1` is the natural choice, but overlapping ranges make the two-pointer merge more complex than for substitutions
- **The `seq` field on child-edge indels can be wrong.** When the Fitch pass creates an indel for an edge B-C, it records the sequence content visible at B. But if positions in that range were already gapped by an earlier edge A-B, the recorded `seq` reflects the gapped state (gaps or stale content), not the original content. Composition must use the parent edge's `seq` for the overlapping portion (cases 1, 6)
- **Insertion-deletion interaction at the same range produces a pair, not a single indel.** The `InDel` type has no "replacement" variant. Case 5 (del + ins, different content) must emit two indels in a defined order. This means the composed output may contain adjacent del+ins pairs at the same range, which consumers (`reconstruct_map_seq`, `add_indel`) must handle in order
- **Composition changes the `indels` field contract.** Currently, indels on an edge are "raw Fitch output" and consumers apply them in list order. After composition, indels become "net annotations" with no redundancy. Verify that all consumers (sequence reconstruction, composition counting, indel count for branch length) produce the same results with the composed list
- **Dense partition indels.** `struct DenseEdgePartition` ([packages/treetime/src/representation/payload/dense.rs#L33](../../packages/treetime/src/representation/payload/dense.rs#L33)) also has an `indels: Vec<InDel>` field at [line 37](../../packages/treetime/src/representation/payload/dense.rs#L37). Verify whether dense partitions undergo edge merge operations that concatenate indels. If so, the same composition logic applies

## Research directions and follow-up

Before implementation:

- Run `flu/h3n2/200` and `ebola` through the `optimize` command with logging at the two concatenation sites. Count: (a) how many indel pairs are non-overlapping (case 7), (b) how many are adjacent (case 2), (c) how many overlap (cases 1, 3-6). This determines whether a minimal implementation (cases 1, 2, 7 only) is sufficient
- Check whether indels on the same edge can have overlapping ranges in the Fitch output. If not, overlapping ranges can only arise from concatenation across edges, simplifying the preconditions

After implementation:

- Property test: `compose_indels(parent, child)` applied to root sequence should produce the same result as applying parent indels then child indels sequentially. This is the round-trip invariant
- Golden master: compare `Composition` counts before and after the change on real datasets. The counts should change (fixing double-counting), so golden masters need updating with justification
- Investigate whether `merge_sibling_pair()` should compose shared indels (analogous to shared substitutions). Currently indels are not shared between siblings, which may miss optimization opportunities in indel-heavy datasets

## Suggested implementation

1. Implement `fn compose_indels(parent: &[InDel], child: &[InDel]) -> Vec<InDel>` in [packages/treetime/src/seq/indel.rs](../../packages/treetime/src/seq/indel.rs). Both input slices should be sorted by `range.0` (add sorting at the call site or as a precondition, matching how `fn compose_substitutions()` requires sorted input). Use a two-pointer merge on `range.0`
2. Add `fn chain_fitch_indels(&self, suffix: &[InDel]) -> Vec<InDel>` to `SparseEdgePartition` ([packages/treetime/src/representation/payload/sparse.rs#L104](../../packages/treetime/src/representation/payload/sparse.rs#L104)), analogous to `fn chain_fitch_subs()` ([packages/treetime/src/representation/payload/sparse.rs#L141](../../packages/treetime/src/representation/payload/sparse.rs#L141))
3. Replace the two concatenation sites (`fn collapse_edge()`, `fn apply_reroot()`) with composition calls
4. Add unit tests covering cases 1-7 above, following the pattern in [packages/treetime/src/seq/\_\_tests\_\_/test_mutation.rs](../../packages/treetime/src/seq/__tests__/test_mutation.rs)

Before implementing cases 3-6, running a dataset (e.g. `flu/h3n2/200`) through the `optimize` command (which performs repeated collapse/reroot cycles) and logging overlap patterns would clarify which cases occur in practice vs remain theoretical.

## Glossary

1. <a id="gloss-1"></a> **CTMC (continuous-time Markov chain).** A stochastic process on a discrete state space where transitions occur in continuous time according to exponential holding times. The transition probability matrix is $P(t) = e^{Qt}$ where $Q$ is the rate matrix. All phylogenetic substitution models (JC69, HKY, GTR) are CTMCs on the nucleotide or amino acid alphabet. [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Semigroup property.** The composition law $P(t+s) = P(t) \cdot P(s)$ satisfied by CTMC transition matrices. Guarantees that the net effect of two consecutive time intervals equals the effect of the combined interval. This property makes substitution composition on consecutive phylogenetic edges exact. [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **TKF91.** The first continuous-time Markov model of sequence evolution incorporating insertions and deletions, introduced by <a id="cite-1b"></a>[Thorne et al. 1991](https://doi.org/10.1007/BF02193625) [[1](#ref-1)]. Restricts indels to single residues, reducing the process to a linear birth-death chain with an exact closed-form pair HMM solution. TKF92 (<a id="cite-2b"></a>[Thorne et al. 1992](https://doi.org/10.1007/BF00163848) [[2](#ref-2)]) extends it to multi-residue indels via indivisible sequence fragments. [↩](#gloss-use-3)
4. <a id="gloss-4"></a> **PIP (Poisson Indel Process).** An indel model where insertion events follow a Poisson process on the phylogeny, reducing the computational complexity of marginal likelihood from exponential (TKF91) to linear in the number of taxa (<a id="cite-3b"></a>[Bouchard-Cote and Jordan 2013](https://www.pnas.org/doi/full/10.1073/pnas.1220450110) [[3](#ref-3)]). [↩](#gloss-use-4)

## References

1. <a id="ref-1"></a> Thorne, Jeffrey L., Hirohisa Kishino, and Joseph Felsenstein. 1991. "An Evolutionary Model for Maximum Likelihood Alignment of DNA Sequences." _Journal of Molecular Evolution_ 33:114-124. https://doi.org/10.1007/BF02193625 [↩¹](#cite-1) [↩²](#cite-1b)
2. <a id="ref-2"></a> Thorne, Jeffrey L., Hirohisa Kishino, and Joseph Felsenstein. 1992. "Inching toward Reality: An Improved Likelihood Model of Sequence Evolution." _Journal of Molecular Evolution_ 34:3-16. https://doi.org/10.1007/BF00163848 [↩¹](#cite-2) [↩²](#cite-2b)
3. <a id="ref-3"></a> Bouchard-Cote, Alexandre, and Michael I. Jordan. 2013. "Evolutionary Inference via the Poisson Indel Process." _Proceedings of the National Academy of Sciences_ 110:1160-1166. https://www.pnas.org/doi/full/10.1073/pnas.1220450110. DOI: [10.1073/pnas.1220450110](https://doi.org/10.1073/pnas.1220450110) [↩¹](#cite-3) [↩²](#cite-3b)
4. <a id="ref-4"></a> Holmes, Ian. 2020. "A Model of Indel Evolution by Finite-State, Continuous-Time Machines." _Genetics_ 216:1187-1204. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7768254/. DOI: [10.1534/genetics.120.303630](https://doi.org/10.1534/genetics.120.303630) [↩](#cite-4)
5. <a id="ref-5"></a> Redelings, Benjamin D., Ian Holmes, Gerton Lunter, and Ari Loytynoja. 2024. "Insertions and Deletions: Computational Methods, Evolutionary Dynamics, and Biological Applications." _Molecular Biology and Evolution_ 41:msae177. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11385596/. DOI: [10.1093/molbev/msae177](https://doi.org/10.1093/molbev/msae177) [↩](#cite-5)
