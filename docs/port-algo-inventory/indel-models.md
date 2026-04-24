# Indel Models

[Back to index](_index.md)

For the full survey of indel modeling approaches in phylogenetics (12 models, 17 software tools, 30 references), see the [Indel Models in Phylogenetics](../reports/indel-models/_index.md) report.

This page covers only the algorithm implemented in v1.

## Poisson Indel Contribution

Adds a Poisson indel log-likelihood term to per-edge branch length optimization. For $k$ observed indel events on a branch of length $t$ with global rate $\mu$: $\ell(t) = k \ln(\mu t) - \mu t - \ln(k!)$. Derivatives $k/t - \mu$ and $-k/t^2$ enter the Newton step alongside substitution derivatives. The rate $\hat{\mu} = \sum_e k_e / \sum_e t_e$ is estimated from the tree at each optimization round.

v1: [`packages/treetime/src/commands/optimize/optimize_indel.rs`](../../packages/treetime/src/commands/optimize/optimize_indel.rs).

v0: not implemented. v0 ignores indels in the likelihood, same as RAxML, IQ-TREE, PhyML.

This is a v1-only feature. See [indel models report](../reports/indel-models/_index.md) for the full catalog of alternative approaches, [intentional change](../port-intentional-changes/optimize-indel-contribution-to-likelihood.md), and [alternatives proposal](../port-proposals/optimize-indel-model-alternatives.md).

## Indel Composition on Consecutive Edges

Substitution models are <a id="gloss-use-1"></a>continuous-time Markov chains <sup>[1](#gloss-1)</sup> (CTMCs) on a finite state space (the alphabet). The transition probability matrix $P(t) = e^{Qt}$ satisfies the <a id="gloss-use-2"></a>semigroup property <sup>[2](#gloss-2)</sup> $P(t+s) = P(t) \cdot P(s)$. This is what makes substitution composition on consecutive edges exact: the net effect of two edges equals the combined transition. `fn compose_substitutions()` ([packages/treetime/src/seq/mutation.rs#L72](../../packages/treetime/src/seq/mutation.rs#L72)) implements this with chaining (A->G + G->T = A->T) and cancellation (A->G + G->A = no change).

The semigroup property does not extend to indel processes. The state space for indel models is the set of all possible sequences (infinite and variable-length), making the matrix exponential intractable. The <a id="gloss-use-3"></a>TKF91 <sup>[3](#gloss-3)</sup> model (<a id="cite-1"></a>[Thorne et al. 1991](https://doi.org/10.1007/BF02193625) [[1](#ref-1)]) is the only indel model with an exact closed-form solution, restricted to single-residue indels. Multi-residue indel models (TKF92 <a id="cite-2"></a>[Thorne et al. 1992](https://doi.org/10.1007/BF00163848) [[2](#ref-2)], the general geometric indel model, the <a id="gloss-use-4"></a>Poisson Indel Process <sup>[4](#gloss-4)</sup> <a id="cite-3"></a>[Bouchard-Cote and Jordan 2013](https://www.pnas.org/doi/full/10.1073/pnas.1220450110) [[3](#ref-3)]) require approximations: pair HMMs, transducer composition, or ODE-based approaches (<a id="cite-4"></a>[Holmes 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7768254/) [[4](#ref-4)]; <a id="cite-5"></a>[Redelings et al. 2024](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11385596/) [[5](#ref-5)]).

The key difference from substitutions: indel composition is not a simple position-wise operation. Insertions change sequence length, so downstream positions shift. Deletions can overlap with prior insertions or other deletions in ways that have no analogue in substitution composition (where each position is independent).

For treetime's purposes, the indel representation is simpler than the general theoretical case because indels are stored as alignment-coordinate ranges (positions do not shift), the Fitch parsimony pass already resolved which positions are gapped vs present, and composition only needs to merge the observed indel annotations, not compute transition probabilities. Treetime's indel composition is a bookkeeping operation on alignment coordinates, not a probabilistic inference step.

Known issue tracking implementation: [M-representation-indel-composition-missing.md](../port-known-issues/M-representation-indel-composition-missing.md).

### Data model

`struct InDel` ([packages/treetime/src/seq/indel.rs#L6](../../packages/treetime/src/seq/indel.rs#L6)):

- `range: (usize, usize)` - half-open position range `[start, end)` in the alignment
- `seq: Seq` - the nucleotide content occupying that range (length = `end - start`)
- `deletion: bool` - `true`: parent has `seq` at this range, child has gaps. `false` (insertion): parent has gaps, child has `seq`

### Composition cases

For two indels on consecutive edges (parent edge first, child edge second):

**Case 1: overlapping deletions (same range)**

```
node A:    ACGTAC    root
node B:    A--TAC    edge A-B: del (1,3) seq="CG"
node C:    A---AC    edge B-C: del (1,4) seq="CGT"  <-- seq is WRONG
expected:  del (1,4) seq="CGT"    from node A (actual pre-deletion content)
```

**Case 2: adjacent deletions (abutting)**

```
node A:    ACGTCG    root
node B:    A--TCG    edge A-B: del (1,3) seq="CG"
node C:    A----G    edge B-C: del (3,5) seq="TC"
expected:  del (1,5) seq="CGTC"
```

**Case 3: overlapping insertions (same range)**

```
node A:    A---AC    root (gaps at 1..4)
node B:    ATTTAC    edge A-B: ins (1,4) seq="TTT"
node C:    AGGGAC    edge B-C: ins (1,4) seq="GGG"
expected:  ins (1,4) seq="GGG"    child sequence wins
```

**Case 4: deletion then insertion (same range, cancellation)**

```
node A:    ACGTAC    root
node B:    A---AC    edge A-B: del (1,4) seq="CGT"
node C:    ACGTAC    edge B-C: ins (1,4) seq="CGT"
expected:  (empty)                indels cancel
```

**Case 5: deletion then insertion (same range, different content)**

```
node A:    ACGTAC    root
node B:    A---AC    edge A-B: del (1,4) seq="CGT"
node C:    ATTTAC    edge B-C: ins (1,4) seq="TTT"
expected:  del (1,4) seq="CGT"  + ins (1,4) seq="TTT"
```

**Case 6: partially overlapping deletions**

```
node A:    ACGTAC    root
node B:    A--TAC    edge A-B: del (1,3) seq="CG"
node C:    A----C    edge B-C: del (2,5) seq="TAC"  <-- seq is WRONG
expected:  del (1,5) seq="CGTA"   merged from parent original content
```

**Case 7: non-overlapping (no interaction)**

```
node A:    ACGTACGT    root
node B:    A--TACGT    edge A-B: del (1,3) seq="CG"
node C:    A--TAC-T    edge B-C: del (6,7) seq="G"
expected:  del (1,3) seq="CG"  + del (6,7) seq="G"
```

### Composition rules summary

| Case | Parent edge    | Child edge     | Result                                        |
| ---- | -------------- | -------------- | --------------------------------------------- |
| 1    | del at range R | del at range R | one del at R with parent's seq                |
| 2    | del + del      | (adjacent)     | merge into one del spanning both ranges       |
| 3    | ins at range R | ins at range R | one ins at R with child's seq                 |
| 4    | del at range R | ins at R, same | cancel (no net indel)                         |
| 5    | del at range R | ins at R, diff | keep both (del + ins, no replacement variant) |
| 6    | del at range R | del at range S | decompose overlap, merge sub-ranges           |
| 7    | any at range R | any at range S | no interaction, both pass through             |

Cases 1, 2, and 7 are expected to dominate in phylogenetic data. Cases 4 and 5 require an indel event to be immediately undone on the next branch, which is biologically unusual. Case 6 requires two indel events at similar but not identical positions on adjacent branches.

### Implementation target

- `fn compose_indels(parent: &[InDel], child: &[InDel]) -> Vec<InDel>` in [packages/treetime/src/seq/indel.rs](../../packages/treetime/src/seq/indel.rs), following the merge-sort pattern of `fn compose_substitutions()`
- `fn chain_fitch_indels()` on `SparseEdgePartition` ([packages/treetime/src/representation/payload/sparse.rs#L104](../../packages/treetime/src/representation/payload/sparse.rs#L104)), analogous to `fn chain_fitch_subs()` ([packages/treetime/src/representation/payload/sparse.rs#L141](../../packages/treetime/src/representation/payload/sparse.rs#L141))

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
