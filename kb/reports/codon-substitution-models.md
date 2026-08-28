# Codon Substitution Models and Ancestral Reconstruction

## Problem

TreeTime reconstructs ancestral amino acid sequences through two independent paths that disagree:

1. **Nucleotide-first:** reconstruct nucleotides under a 4-state model, then translate the inferred codons.
2. **Amino-acid-first:** translate tip sequences to amino acids, then reconstruct under a 20-state (or 22-state) model.

These yield different ancestral amino acid assignments because the genetic code is a many-to-one map (61 sense codons to 20 amino acids), and composing a nonlinear map with an argmax operator is not commutative. For a map $f$ from codons to amino acids and reconstruction operator $R$:

$$f(R(x)) \neq R(f(x))$$

The augur pipeline's `--report-inconsistent-translation` flag detects the disagreement, confirming that developers are aware of it. The published phylogenetics literature contains no formal analysis of this inconsistency -- instead, the field solved it by reconstructing in codon space, where the problem cannot arise.

## The two paths in the Nextstrain ecosystem

The Nextstrain pipeline (augur + TreeTime + Nextclade) computes ancestral amino acid sequences via two routes that use different operation orders, different models, and different translation implementations.

### Path A: `augur translate` (nucleotide-first)

Reconstruct nucleotides on the tree with TreeTime's `TreeAnc` (`alphabet='nuc'`, JC69, joint mode), store the full nucleotide sequence at each node, then translate each node's CDS with BioPython's `safe_translate` [[src](https://github.com/nextstrain/augur/blob/7c7ac4036a5a5ce06ef2a637d917983ca3b17a40/augur/translate.py#L268-L309)]. The ancestral nucleotide reconstruction is triggered via `augur ancestral` [[src](https://github.com/nextstrain/augur/blob/7c7ac4036a5a5ce06ef2a637d917983ca3b17a40/augur/ancestral.py#L272-L321)].

### Path B: `augur ancestral --translations` (amino-acid-first, used by flu)

Translate tip nucleotide sequences to amino acids with Nextclade's codon-aware translator (reference-coordinate peptides), then run a second, independent `TreeAnc` reconstruction directly in amino acid space (`alphabet='aa'`, JC69 with 22 states, joint mode) [[src](https://github.com/nextstrain/augur/blob/7c7ac4036a5a5ce06ef2a637d917983ca3b17a40/augur/ancestral.py#L580-L685)]. Augur's `_validate_translated_consistency` [[src](https://github.com/nextstrain/augur/blob/7c7ac4036a5a5ce06ef2a637d917983ca3b17a40/augur/ancestral.py#L523-L567)] detects disagreements between the two paths.

### Four root causes of disagreement

**RC1 -- Non-commutativity.** Reconstructing three nucleotide sites independently then translating is not the same as reconstructing one amino acid. The nucleotide path is blind to codon structure: it can resolve a per-site tie in a way that produces a codon encoding an amino acid present in no child node.

**RC2 -- Different models.** Path A uses nucleotide JC69 (4-state). Path B uses amino acid JC69 (uniform 22-state, including gap and stop as substitutable states) [[src](../../packages/legacy/treetime/treetime/nuc_models.py#L18-L40)]. Neither is a codon model; both are approximations biased in different directions.

**RC3 -- Two different translators.** Path A uses BioPython's `safe_translate`; Path B uses Nextclade's codon-aware translator. They diverge on the hard cases: frameshifts, indels, and ambiguous codons.

**RC4 -- Encoding defects in the AA model.** TreeTime's amino acid alphabet treats `-` (gap) and `*` (stop) as substitutable states [[src](../../packages/legacy/treetime/treetime/seq_utils.py#L22)], and defines `X` as all-ones over all 22 symbols including gap and stop [[src](../../packages/legacy/treetime/treetime/seq_utils.py#L92)].

### Worked examples

**E1 -- Synonymous tie-break creates a phantom residue.** Three leaves under an internal node, codon at one position:

```
leaf1  CTG = Leu
leaf2  TTG = Leu
leaf3  ATT = Ile
```

Path B (AA space): internal states {Leu, Leu, Ile} -> reconstructs **Leu** (majority). Path A (nuc space): position 1 has three distinct bases {C, T, A}, each appearing once -- a three-way JC69 tie. If the tie-break yields `A`, the internal codon becomes `ATG` = **Met**, an amino acid present in no child. The nucleotide reconstruction is correct per its model but wrong at the protein level.

**E2 -- Frameshift window.** A tip with a single-base deletion causing a 30-nt frameshift. Nextclade detects the frameshift, masks affected nucleotides to `N`, fills the peptide window with `X` [[src](https://github.com/nextstrain/nextclade/blob/629ea8fb0cfe3165db13ddc75dab16d3c2292dbd/packages/nextclade/src/translate/translate_genes.rs#L191-L224)]. `safe_translate` translates whatever codons the annotation extracts at the shifted frame offset -- producing real-looking but incorrect amino acids. Path A feeds those spurious residues into reconstruction as if they were signal.

**E3 -- Ambiguous codon `ACN`.** All four resolutions (ACA, ACC, ACG, ACT) encode threonine. Nextclade enumerates -> emits `T` (correct, informative) [[src](https://github.com/nextstrain/nextclade/blob/629ea8fb0cfe3165db13ddc75dab16d3c2292dbd/packages/nextclade/src/translate/translate.rs#L25)]. `safe_translate` fails the table lookup for `ACN` -> emits `X` (discards recoverable information). Here Path A is the weaker translator.

### Defect catalog

**Algorithmic:**

- **D1** -- No codon model anywhere. Both paths ignore codon structure and synonymous/nonsynonymous asymmetry.
- **D2** -- Path B uses uniform 22-state JC69 for amino acids. Real AA substitution is non-uniform (LG/WAG/BLOSUM).

**TreeTime AA alphabet encoding:**

- **D3** -- `-` (gap) and `*` (stop) are substitutable states in the AA rate matrix. An internal node can be reconstructed as "an evolving gap."
- **D4** -- `X` is all-ones over all 22 symbols including gap and stop. Correct `X` should span only the 20 real amino acids.
- **D5** -- `seq2array` defaults `ambiguous='N'`; `N` is not an AA symbol, causing silent misencoding. The v1 Rust port must use `X` as the AA ambiguous character.

**Translator disconnect:**

- **D6** -- Two translators (BioPython vs Nextclade) with divergent gap, frameshift, insertion, and ambiguity semantics.
- **D7** -- Nextclade strips insertions relative to the reference [[src](https://github.com/nextstrain/nextclade/blob/629ea8fb0cfe3165db13ddc75dab16d3c2292dbd/packages/nextclade/src/align/insertions_strip.rs#L54-L103)]. Path B cannot see amino acid changes carried by in-frame insertions.
- **D8** -- Frameshift body and genuine protein-level ambiguity both map to `X` in Nextclade. Downstream reconstruction cannot distinguish "misaligned garbage" from "uncertain residue."
- **D9** -- Reference/root mismatch risk: Path A derives root AA from `translate(reference.nuc)`; Path B may use `--aa-root-sequence`, `translate(nuc ref)`, or an inferred root. If the Nextclade dataset reference differs from the annotation reference, mutation polarity diverges.

**Robustness:**

- **D10** -- Length coupling enforced by `assert`, not reconciled. Annotation/CDS-length mismatch causes a crash.
- **D11** -- `safe_translate` errors when CDS length is not divisible by 3; Nextclade silently floors to a whole codon.

## Codon substitution models

Codon models describe evolution of protein-coding sequences using the codon as the unit of evolution. The state space is the 61 sense codons under the standard genetic code; the three stop codons are excluded. The instantaneous rate matrix $Q$ is $61 \times 61$.

Two foundational models were proposed independently in 1994:

### MG94

<a id="cite-1"></a>[Muse and Gaut 1994](https://doi.org/10.1093/oxfordjournals.molbev.a040152) [[1](#ref-1)]. Substitution rate between codons $i$ and $j$ differing at a single nucleotide position:

$$q_{ij} = \begin{cases} \mu \, \pi_{n(j)} & \text{synonymous} \\ \mu \, \omega \, \pi_{n(j)} & \text{nonsynonymous} \end{cases}$$

where:

- $\pi_{n(j)}$ -- equilibrium frequency of the target nucleotide (the changed position in codon $j$)
- $\omega$ -- nonsynonymous/synonymous rate ratio ($d_N/d_S$)
- $\mu$ -- overall mutation rate

Codons differing by more than one nucleotide have $q_{ij} = 0$.

The MG parameterization models mutation at the nucleotide level. Mutation-selection models (FMutSel, <a id="cite-6"></a>[Yang and Nielsen 2008](https://doi.org/10.1093/molbev/msm284) [[6](#ref-6)]) were built on this formulation.

### GY94

<a id="cite-2"></a>[Goldman and Yang 1994](https://doi.org/10.1093/oxfordjournals.molbev.a040153) [[2](#ref-2)]. Substitution rate between codons $i$ and $j$ differing at a single nucleotide position:

$$q_{ij} = \begin{cases} \pi_j \, \kappa & \text{synonymous transition} \\ \pi_j & \text{synonymous transversion} \\ \pi_j \, \kappa \, \omega & \text{nonsynonymous transition} \\ \pi_j \, \omega & \text{nonsynonymous transversion} \end{cases}$$

where:

- $\pi_j$ -- equilibrium frequency of the target codon
- $\kappa$ -- transition/transversion rate ratio
- $\omega$ -- nonsynonymous/synonymous rate ratio

GY parameterizes rates proportional to the target codon frequency rather than the target nucleotide frequency. This is less mechanistic but more flexible for protein-level constraints. GY also incorporates the transition/transversion bias $\kappa$ and optionally a physicochemical distance between the encoded amino acids.

### MG vs GY

The key distinction: MG uses target nucleotide frequency $\pi_{n(j)}$; GY uses target codon frequency $\pi_j$. Bayesian comparison <a id="cite-7"></a>[Rodrigue et al. 2008](https://doi.org/10.1534/genetics.108.092254) [[7](#ref-7)] favored MG-style parameterization: it maintains distinct nucleotide-level mutation parameters while allowing flexible amino acid or codon preferences.

Both produce a $61 \times 61$ rate matrix $Q$, normalized so that the expected rate equals 1:

$$\beta = -\sum_i \pi_i \, Q_{ii}$$

The transition probability matrix over branch length $t$ is $P(t) = e^{Qt/\beta}$.

### Extensions

**Site-heterogeneous models.** <a id="cite-3"></a>[Yang et al. 2000](https://doi.org/10.1093/genetics/155.1.431) [[3](#ref-3)] introduced models allowing $\omega$ to vary among amino acid sites (M0-M8 in PAML's `codeml`). This enables detection of positive selection at specific sites even when gene-wide $\omega < 1$.

**Empirical codon models.** <a id="cite-4"></a>[Schneider et al. 2005](https://doi.org/10.1186/1471-2105-6-134) [[4](#ref-4)] estimated the first empirical $61 \times 61$ matrix from vertebrate genomes. <a id="cite-5a"></a>[Kosiol et al. 2007](https://doi.org/10.1093/molbev/msm064) [[5](#ref-5)] built the ECM (Empirical Codon Model), which incorporates instantaneous multi-nucleotide changes (doublet and triplet substitutions) and consistently outperforms mechanistic models. <a id="cite-8"></a>[Doron-Faigenboim and Pupko 2007](https://doi.org/10.1093/molbev/msl175) [[8](#ref-8)] combined empirical amino acid replacement rates with mechanistic codon parameters ($\omega$, $\kappa$).

**Mutation-selection models.** FMutSel [[6](#ref-6)] separates mutation bias from selection on codon usage. Codon fitness parameters together with mutation-bias parameters predict optimal codon frequencies. Population-genetic grounding allows estimation of selective strengths.

## Ancestral reconstruction under codon models

Ancestral reconstruction under codon models uses the same Felsenstein pruning algorithm described in [ancestral.md](../algo/ancestral.md), with the state space set to 61 codons instead of 4 nucleotides or 20 amino acids.

### Backward pass (leaf to root)

For internal node $k$ with children $i$ and $j$, the partial likelihood for codon $c$ at site $s$:

$$w_k(c) = \left[\sum_{c'} P(c \to c' \mid t_i) \, w_i(c')\right] \cdot \left[\sum_{c''} P(c \to c'' \mid t_j) \, w_j(c'')\right]$$

where $P(c \to c' \mid t)$ is the $(c, c')$ entry of the $61 \times 61$ transition probability matrix $e^{Qt}$.

### Forward pass (root to leaf)

Identical to the nucleotide case: each node receives an outgroup message via cavity/division, combined with the backward message to produce the marginal posterior over 61 codons at each site.

### Deriving amino acid and nucleotide sequences

After codon reconstruction, amino acid and nucleotide sequences are deterministic projections:

**Amino acid posterior:** For amino acid $a$, sum codon posteriors over all synonymous codons:

$$P(a \mid \text{data}) = \sum_{c : \text{translate}(c) = a} P(c \mid \text{data})$$

**Nucleotide triplet:** The MAP codon directly encodes the three nucleotides.

This derivation produces consistent nucleotide and amino acid ancestral sequences from a single reconstruction. The non-commutativity problem cannot arise because there is no separate nucleotide or amino acid reconstruction step.

## How established tools implement codon reconstruction

Source inspection of three major phylogenetic tools confirms a convergent architecture.

### PAML/codeml

`codeml.c` [[src](https://github.com/abacus-gene/paml/blob/4c7902fe972737ef5e80bb18f159a6e6acace3d3/src/codeml.c)], `treesub.c` [[src](https://github.com/abacus-gene/paml/blob/4c7902fe972737ef5e80bb18f159a6e6acace3d3/src/treesub.c)].

State space: `com.ncode = Nsensecodon` (61). Both marginal (`AncestralMarginal()`) and joint (`AncestralJointPPSG2000()`, <a id="cite-9"></a>[Pupko et al. 2000](https://doi.org/10.1093/oxfordjournals.molbev.a026369) [[9](#ref-9)]) reconstruction operate on the 61-codon state space. Amino acid probabilities are derived by summing codon posteriors over synonymous codons. Stop codons are excluded via `FROM61[]`/`FROM64[]` lookup tables.

### IQ-TREE 2

`modelcodonparametric.cpp` [[src](https://github.com/iqtree/iqtree2/blob/a00094e03d1ae984e1497e16738f91514df8c366/model/modelcodonparametric.cpp)], `phylotreesse.cpp` [[src](https://github.com/iqtree/iqtree2/blob/a00094e03d1ae984e1497e16738f91514df8c366/tree/phylotreesse.cpp)].

Implements MG, GY, ECM, and ECMrest codon models <a id="cite-10"></a>[Minh et al. 2020](https://doi.org/10.1093/molbev/msaa015) [[10](#ref-10)]. `num_states = 61` when a codon model is active. Ancestral reconstruction (`computeMarginalAncestralState`) is model-agnostic -- the same code path handles 4-state (DNA), 20-state (AA), and 61-state (codon) models. GY is the default codon model. Joint reconstruction (Pupko et al. 2000) is implemented but currently disabled; only marginal is active.

### HyPhy

`likefunc.cpp` [[src](https://github.com/veg/hyphy/blob/646e16eb4243e5e5588ac340483fcb987552751d/src/core/likefunc.cpp)], `ancestral.bf` [[src](https://github.com/veg/hyphy/blob/646e16eb4243e5e5588ac340483fcb987552751d/res/TemplateBatchFiles/libv3/tasks/ancestral.bf)].

`ReconstructAncestors()` and `SampleAncestors()` operate on the state space of the active likelihood function. With a codon model, `CHARS = 61` and the ancestral state matrix maps node indices to codon indices. The batch language architecture separates model definition from reconstruction algorithm, but reconstruction is native to the model's state space.

### Common pattern

None of these tools implements "reconstruct nucleotides then translate" or "translate then reconstruct amino acids." The codon model unifies both levels: reconstructed codons yield both nucleotide and amino acid sequences by deterministic projection.

## Computational cost

The rate matrix for a codon model is $61 \times 61$, compared to $4 \times 4$ for nucleotides or $20 \times 20$ for amino acids. Matrix exponentiation via eigendecomposition scales as $O(k^3)$ where $k$ is the number of states. The ratio $(61/4)^3 \approx 3500$ gives the approximate cost increase over nucleotide models per branch.

Total complexity for ancestral reconstruction: $O(n \cdot k^2 \cdot L)$ for $n$ nodes, $k$ states, and $L$ alignment positions (see [ancestral.md](../algo/ancestral.md)). With $k = 61$ this is $(61/4)^2 \approx 230$ times more expensive than nucleotide reconstruction per site.

<a id="cite-5b"></a>[Ren et al. 2005](https://doi.org/10.1080/10635150500354688) [[5](#ref-5)] found codon models "unfeasible for tree search in large data sets." CodonPhyML <a id="cite-11"></a>[Gil et al. 2013](https://doi.org/10.1093/molbev/mst034) [[11](#ref-11)] is the only dedicated fast ML phylogeny tool under codon models.

For TreeTime's use case -- fixed tree topology, variable branch lengths, ancestral reconstruction -- the computational cost is more manageable. Tree search is not needed. The eigendecomposition of $Q$ is computed once and cached; per-branch cost is the matrix-vector product $P(t) \cdot w$, which is $O(k^2)$ per site.

## Information loss from translation

<a id="cite-5c"></a>[Ren et al. 2005](https://doi.org/10.1080/10635150500354688) [[5](#ref-5)] quantified the phylogenetic information lost by translating coding sequences to amino acids before analysis. Using 106 protein-coding genes from 8 yeast species: 30% loss for deep divergences, 66% for recent divergences. Synonymous substitutions (discarded by translation) carry most of the signal at short evolutionary distances and saturate first.

Nucleotide models performed best for recent divergences; amino acid models for deep divergences. Codon models combined the advantages of both, performing well at all evolutionary distances.

## Implications for TreeTime v1

### The current two-path design is an anomaly

No other major phylogenetic tool implements dual reconstruction paths for the same coding data. The field's established solution -- codon-level reconstruction -- predates the augur pipeline and eliminates the inconsistency by construction.

### Codon models are feasible for TreeTime's use case

TreeTime does not perform tree search, which is the primary bottleneck for codon models. The reconstruction algorithm is the same Felsenstein pruning already implemented in v1 ([`packages/treetime/src/partition/marginal_dense.rs`](../../packages/treetime/src/partition/marginal_dense.rs), [`packages/treetime/src/partition/marginal_passes.rs`](../../packages/treetime/src/partition/marginal_passes.rs)). Only the state space and rate matrix change.

### Implementation path

1. Define a 61-state codon alphabet (sense codons, standard genetic code)
2. Implement MG94 or GY94 rate matrix as a $61 \times 61$ $Q$ matrix in [`packages/treetime/src/gtr/`](../../packages/treetime/src/gtr/)
3. Extend the marginal reconstruction passes to operate with `num_states = 61`
4. Derive amino acid posteriors by summing codon posteriors over synonymous codons
5. Derive nucleotide sequences by triplet extraction from MAP codons

The dense implementation is most directly applicable: sparse representation assumes a compact alphabet (4 or 20 states) and Fitch compression, which would need adaptation for 61 states.

### Partition model alternative

Nucleotide models with separate rate parameters for first, second, and third codon positions capture much of the synonymous/nonsynonymous rate asymmetry at lower computational cost. This does not eliminate the non-commutativity problem (nucleotides are still reconstructed independently per position), but provides a cheaper intermediate.

## Open questions

- No published benchmark compares codon-level ancestral state accuracy against nucleotide-then-translate. The Ren et al. 2005 information-loss numbers apply to tree topology, not ancestral states.
- Multi-nucleotide substitutions (doublet/triplet changes) improve model fit [[5](#ref-5)] but add complexity to the $Q$ matrix. Whether they matter for ancestral reconstruction accuracy (as opposed to selection inference) is less studied.
- Codon models assume a fixed alignment and do not model insertions or deletions. Indels relative to the reference are a major source of disagreement between the two translation paths in augur/TreeTime (frameshift handling). Codon models do not solve this; alignment preprocessing must handle it before reconstruction.
- Whether 61-state marginal reconstruction is tractable for TreeTime's largest datasets (thousands of tips, genome-length alignments) needs benchmarking. The sparse representation may need adaptation for 61-state alphabets.

## Glossary

1. <a id="gloss-1"></a> **Sense codon.** A codon that encodes an amino acid. The standard genetic code has 61 sense codons; the remaining 3 (UAA, UAG, UGA) are stop codons. [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Degeneracy.** The property of the genetic code whereby multiple codons encode the same amino acid. The forward map (codon -> amino acid) is many-to-one; the reverse (amino acid -> codon) is one-to-many. [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **ECM (Empirical Codon Model).** A $61 \times 61$ substitution matrix estimated from large-scale alignments of coding sequences, as opposed to mechanistic models that parameterize rates from first principles. [↩](#gloss-use-3)
4. <a id="gloss-4"></a> **FMutSel.** Mutation-selection codon model. Separates mutation bias (nucleotide-level) from selection on codon usage (fitness parameters). Built on the MG formulation. [↩](#gloss-use-4)

## References

1. <a id="ref-1"></a> Muse, Spencer V., and Brandon S. Gaut. 1994. "A Likelihood Approach for Comparing Synonymous and Nonsynonymous Nucleotide Substitution Rates, with Application to the Chloroplast Genome." _Molecular Biology and Evolution_ 11:715-724. https://doi.org/10.1093/oxfordjournals.molbev.a040152 [↩](#cite-1)
2. <a id="ref-2"></a> Goldman, Nick, and Ziheng Yang. 1994. "A Codon-Based Model of Nucleotide Substitution for Protein-Coding DNA Sequences." _Molecular Biology and Evolution_ 11:725-736. https://doi.org/10.1093/oxfordjournals.molbev.a040153 [↩](#cite-2)
3. <a id="ref-3"></a> Yang, Ziheng, Rasmus Nielsen, Nick Goldman, and Anne-Mette Krabbe Pedersen. 2000. "Codon-Substitution Models for Heterogeneous Selection Pressure at Amino Acid Sites." _Genetics_ 155:431-449. https://doi.org/10.1093/genetics/155.1.431 [↩](#cite-3)
4. <a id="ref-4"></a> Schneider, Adrian, Gina M. Cannarozzi, and Gaston H. Gonnet. 2005. "Empirical Codon Substitution Matrix." _BMC Bioinformatics_ 6:134. https://doi.org/10.1186/1471-2105-6-134 [↩](#cite-4)
5. <a id="ref-5"></a> Ren, Fengrong, Hiroshi Tanaka, and Ziheng Yang. 2005. "An Empirical Examination of the Utility of Codon-Substitution Models in Phylogeny Reconstruction." _Systematic Biology_ 54:808-818. https://doi.org/10.1080/10635150500354688 [↩¹](#cite-5a) [↩²](#cite-5b) [↩³](#cite-5c)
6. <a id="ref-6"></a> Yang, Ziheng, and Rasmus Nielsen. 2008. "Mutation-Selection Models of Codon Substitution and Their Use to Estimate Selective Strengths on Codon Usage." _Molecular Biology and Evolution_ 25:568-579. https://doi.org/10.1093/molbev/msm284 [↩](#cite-6)
7. <a id="ref-7"></a> Rodrigue, Nicolas, Nicolas Lartillot, and Herve Philippe. 2008. "Bayesian Comparisons of Codon Substitution Models." _Genetics_ 180:1579-1591. https://doi.org/10.1534/genetics.108.092254 [↩](#cite-7)
8. <a id="ref-8"></a> Doron-Faigenboim, Adi, and Tal Pupko. 2007. "A Combined Empirical and Mechanistic Codon Model." _Molecular Biology and Evolution_ 24:388-397. https://doi.org/10.1093/molbev/msl175 [↩](#cite-8)
9. <a id="ref-9"></a> Pupko, Tal, Itsik Pe'er, Ron Shamir, and Dan Graur. 2000. "A Fast Algorithm for Joint Reconstruction of Ancestral Amino Acid Sequences." _Molecular Biology and Evolution_ 17:890-896. https://doi.org/10.1093/oxfordjournals.molbev.a026369 [↩](#cite-9)
10. <a id="ref-10"></a> Minh, Bui Quang, Heiko A. Schmidt, Olga Chernomor, Dominik Schrempf, Michael D. Woodhams, Arndt von Haeseler, and Robert Lanfear. 2020. "IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era." _Molecular Biology and Evolution_ 37:1530-1534. https://doi.org/10.1093/molbev/msaa015 [↩](#cite-10)
11. <a id="ref-11"></a> Gil, Manuel, Marcelo Serrano Zanetti, Stefan Zoller, and Maria Anisimova. 2013. "CodonPhyML: Fast Maximum Likelihood Phylogeny Estimation under Codon Substitution Models." _Molecular Biology and Evolution_ 30:1270-1280. https://doi.org/10.1093/molbev/mst034 [↩](#cite-11)
