# Indel Models for Phylogenetic Branch Length Optimization

[Back to index](_index.md)

Standard phylogenetic substitution models (JC69, K80, GTR) treat sequence length as fixed and model only character-state changes at individual sites. Insertions and deletions (indels) are either ignored entirely (gaps treated as missing data) or handled by ad hoc gap penalties during alignment, separate from the evolutionary model. This document catalogs approaches to incorporating indel information into the likelihood function, with focus on branch length optimization.

The approaches range from simple heuristics (affine gap penalties) through count-based models (Poisson) to full mechanistic birth-death processes (TKF91, TKF92, PIP). Each makes different tradeoffs between biological realism, computational cost, and implementation complexity.

## Treating gaps as missing data (standard practice)

Most phylogenetic ML software (RAxML, IQ-TREE, PhyML, BEAST, MrBayes, TreeTime v0) treats gap characters as missing data: alignment columns containing gaps contribute no information to the likelihood at the gapped positions. The site likelihood at a gap position marginalizes over all possible states, which is equivalent to ignoring that position.

This is the simplest and most widely used approach. It requires no indel model and adds no parameters. The limitation is that indel events carry phylogenetic signal that is discarded: a branch with zero substitutions but one or more indels is assigned zero length even though the indel represents genuine evolutionary change.

<a id="cite-1"></a>[Warnow 2012](https://doi.org/10.1371/currents.rrn1308) [[1](#ref-1)] proved that ML phylogeny estimation treating indels as missing data is statistically inconsistent, even given the true alignment. This motivates proper evolutionary indel models.

v0: gaps as missing data.
v1: gaps as missing data for substitution likelihood; Poisson indel contribution added (see below).

---

## Affine gap penalty model

Assigns a fixed log-likelihood cost per indel event (gap opening penalty $d$) plus a per-position extension cost ($e$), where $d > e$. Introduced by <a id="cite-2"></a>[Gotoh 1982](<https://doi.org/10.1016/0022-2836(82)90398-9>) [[2](#ref-2)] for alignment scoring; the $O(MN)$ DP algorithm remains a standard alignment technique.

### Mathematical formulation

For a gap of length $g$: $\text{cost}(g) = d + (g - 1) \cdot e$

This implies a geometric distribution on gap lengths, with the ratio $e/d$ controlling the expected length. The total gap penalty for an edge is the sum over all gap events.

### Relationship to evolutionary models

<a id="cite-3"></a>[Rivas 2005](https://doi.org/10.1186/1471-2105-6-63) [[3](#ref-3)] developed time-dependent rate matrix models compatible with affine gap penalties and proved a key limitation: single insertion events with geometric instantaneous length distributions do NOT produce geometric insert lengths at finite evolutionary times. The affine gap model is an approximation without a consistent continuous-time evolutionary process behind it.

### Indel length distributions

Empirical studies show indel lengths follow power-law (Zipf) or multi-exponential distributions, not the geometric distribution implied by affine gaps:

- <a id="cite-4"></a>[Qian and Goldstein 2001](https://doi.org/10.1002/prot.1129) [[4](#ref-4)]: gap lengths from structural alignments are multi-exponential (four components)
- <a id="cite-5"></a>[Cartwright 2009](https://doi.org/10.1093/molbev/msn275) [[5](#ref-5)]: Zipf power-law model fits far better than geometric
- <a id="cite-6"></a>[Wygoda et al. 2024](https://doi.org/10.1093/bioinformatics/btae043) [[6](#ref-6)]: ABC model-selection confirms Zipf fits most empirical datasets better than geometric

### Practical considerations

- Simple to implement: two parameters ($d$, $e$)
- Not probabilistic: does not define a proper likelihood function
- Cannot be combined with substitution likelihood in a principled way
- Widely available in alignment software but not in tree likelihood computation
- The geometric length assumption is empirically wrong

### Software

All major alignment tools (MAFFT, MUSCLE, ClustalW, PRANK) use affine or related gap penalties for alignment, but none incorporate them into the phylogenetic likelihood for branch length estimation.

### Status

v0: not used in likelihood. v1: not used in likelihood.

---

## Poisson indel count model (implemented in v1)

Models the number of indel events on each edge as an independent Poisson process with a single global rate $\mu$. Each indel event contributes equally regardless of length or direction (insertion vs deletion).

### Mathematical formulation

For a branch with $k$ indel events and length $t$:

$$\ell_{\text{indel}}(t) = k \ln(\mu t) - \mu t - \ln(k!)$$

First and second derivatives:

$$\frac{d\ell}{dt} = \frac{k}{t} - \mu, \qquad \frac{d^2\ell}{dt^2} = -\frac{k}{t^2}$$

The global rate estimate is the Poisson MLE:

$$\hat{\mu} = \frac{\sum_e k_e}{\sum_e t_e}$$

where the sums run over all edges $e$ in the tree.

### Key properties

- Strictly concave in $t$ for $k > 0$ (guaranteed unique optimum)
- At $t = 0$ with $k > 0$: $\ell \to -\infty$ and $d\ell/dt \to +\infty$ (forces positive branch length)
- MLE for branch length from indels alone: $t^* = k / \mu$
- The indel contribution is additive with the substitution log-likelihood in the Newton step

### Assumptions and limitations

- Each indel event has equal weight regardless of indel length (a 100-base deletion counts the same as a 1-base deletion)
- Insertions and deletions are treated identically (single rate, not separate birth-death processes)
- The rate is constant across the tree (no branch-specific variation)
- No mechanistic model of the indel process; does not distinguish insertion from deletion or model indel lengths
- Not intended for reconstructing the indel process or estimating indel-specific evolutionary parameters

### Practical considerations

- Minimal implementation effort: one function plus wiring into existing Newton optimizer
- No new user-facing parameters (rate estimated automatically)
- Computational cost negligible: one scalar per edge per iteration
- Integrates into eigendecomposition-based branch length optimization without changing the optimization structure
- Appropriate when the goal is preventing zero-length assignment on indel-only branches, not reconstructing indel history

### Software

No other phylogenetic software uses this approach. It is a TreeTime v1 design-doc feature. Most software ignores indels entirely in the likelihood.

### Status

v0: not implemented.
v1: implemented in [`packages/treetime/src/commands/optimize/optimize_indel.rs`](../../packages/treetime/src/commands/optimize/optimize_indel.rs). See [intentional change](../port-intentional-changes/optimize-indel-contribution-to-likelihood.md) and [design doc](../algorithms/optimize.md#poisson-indel-contribution-implemented).

### Related known issues

- [Grid search zero-comparison ignores indel likelihood](../port-known-issues/M-optimize-grid-zero-ignores-indels.md)
- [Timetree branch length distribution ignores indels](../port-known-issues/N-timetree-branch-length-distribution-ignores-indels.md)

---

## TKF91 (Thorne-Kishino-Felsenstein 1991)

The first continuous-time Markov chain model that jointly models substitutions and indels as a single evolutionary process. Introduced by <a id="cite-7"></a>[Thorne, Kishino, and Felsenstein 1991](https://doi.org/10.1007/BF02193625) [[7](#ref-7)].

### Mathematical formulation

Sequence evolution is modeled as a birth-death process on a linear chain of characters:

- Each existing character can spawn a new adjacent character (insertion) at rate $\lambda$
- Each existing character can be removed (deletion) at rate $\mu$
- Constraint: $\lambda < \mu$ for stationarity
- Substitutions at each site follow a standard model (JC69, HKY85, GTR, etc.)
- An "immortal link" anchors the left end of the sequence (prevents the entire sequence from being deleted)

The equilibrium sequence length distribution is geometric with parameter $1 - \lambda/\mu$.

Key transition probabilities for branch length $t$:

$$\beta(t) = \frac{1 - e^{(\lambda - \mu)t}}{\mu - \lambda e^{(\lambda - \mu)t}}$$

$$p_1(t) = e^{-\mu t}(1 - \lambda \beta) \quad \text{(homolog survival)}$$

$$p_0'(t) = \mu \beta \quad \text{(deletion probability)}$$

The model can be represented as a three-state pair HMM (match, insert, delete) with transition probabilities derived from $\lambda$, $\mu$, and $t$.

### Computational complexity

- Pairwise alignment: $O(L_1 \cdot L_2)$ via pair HMM DP
- Multiple sequences on a tree: $O(L^N)$ exact (exponential in number of taxa $N$)
- With fixed alignment: $O(L \cdot N)$ postorder traversal (Felsenstein peeling extended to handle indel events per column)

### Key limitations

- Single-residue indels only (each insertion/deletion affects exactly one character). Biologically unrealistic: real indels commonly span multiple positions.
- Site independence within the birth-death process
- Exponential multi-sequence complexity without approximations

### Practical considerations for branch length optimization

- Two parameters ($\lambda$, $\mu$) vs one for Poisson ($\mu$)
- Distinguishes insertions from deletions (separate rates)
- Requires modifying the tree likelihood computation structure (pair HMM or extended peeling)
- Cannot be simply added to the existing Newton step; requires restructuring the per-edge likelihood evaluation
- Over-parameterized for the branch-length-prevents-zero use case

### Software implementations

| Software   | Language | Purpose                            | Reference                              |
| :--------- | :------- | :--------------------------------- | :------------------------------------- |
| BAli-Phy   | C++      | Bayesian joint alignment/phylogeny | Redelings and Suchard 2005             |
| BEAST      | Java     | Bayesian MCMC phylogenetics        | Drummond et al. (beast-dev/beast-mcmc) |
| StatAlign  | Java     | Statistical alignment with MCMC    | Novak et al. 2008                      |
| rust-phylo | Rust     | Tree likelihood                    | acg-team/rust-phylo                    |
| CRAN/TKF   | C/R      | Pairwise distance estimation       | Tan (CRAN package)                     |

### Status

v0: not implemented. v1: not implemented.

### References

- <a id="cite-7-loc"></a>[Thorne, Kishino, and Felsenstein 1991](https://doi.org/10.1007/BF02193625) (TKF91 original paper)
- <a id="cite-8"></a>[Holmes 2005](https://doi.org/10.1093/bioinformatics/bti177) [[8](#ref-8)] (EM estimation of TKF91 insertion/deletion rates)
- <a id="cite-9"></a>[Rivas 2008](https://doi.org/10.1371/journal.pcbi.1000172) [[9](#ref-9)] (extended Felsenstein peeling for single-residue indels)

---

## TKF92 (Thorne-Kishino-Felsenstein 1992)

Extension of TKF91 to multi-residue indels. Introduced by <a id="cite-10"></a>[Thorne, Kishino, and Felsenstein 1992](https://doi.org/10.1007/BF00163848) [[10](#ref-10)].

### Mathematical formulation

Sequences are composed of indivisible "fragments" (blocks of linked residues). Insertions and deletions operate on entire fragments. Fragment lengths follow a geometric distribution. Parameters: $\lambda$ (insertion rate), $\mu$ (deletion rate), plus a fragment length parameter.

TKF91 is the special case where all fragments have length 1.

### Computational complexity

Pairwise: $O(L_1 \cdot L_2)$. Multi-sequence: approximate (no exact polynomial solution for arbitrary trees).

### Key difference from TKF91

Allows multi-residue indels, addressing the main biological limitation of TKF91. The geometric fragment length distribution is more realistic than single-residue-only but still does not match the empirical power-law distribution.

### Status

v0: not implemented. v1: not implemented.

---

## Long indel model (SID model, Miklos-Lunter-Holmes 2004)

A continuous-time Markov model allowing indels of arbitrary length with geometric length distributions. Introduced by <a id="cite-11"></a>[Miklos, Lunter, and Holmes 2004](https://doi.org/10.1093/molbev/msh043) [[11](#ref-11)].

### Mathematical formulation

Extends TKF91/92 by allowing single insertion and deletion events of any length, with geometric length distributions (parameterized separately for insertions and deletions). Introduces "rate grammar" notation and uses an infinite sequence embedding to make rates position-independent.

The key algorithmic contribution is the "chop zone" decomposition: the alignment is partitioned into zones where the trajectory likelihood can be computed independently, yielding $O(N^2)$ pairwise complexity.

### Parameters

Insertion rate $\lambda$, deletion rate $\mu$, mean insertion length $\bar{X}$, mean deletion length $\bar{Y}$. TKF91 is the special case $\bar{X} = \bar{Y} = 1$.

### Practical considerations

- More biologically realistic than TKF91/92 (arbitrary-length indels)
- No exact closed-form solution for finite-time transition probabilities; requires numerical ODE integration or finite-state machine approximation
- <a id="cite-12"></a>[Holmes 2020](https://doi.org/10.1534/genetics.120.303630) [[12](#ref-12)] provides systematic approximation via coarse-graining of pair HMM state spaces
- Too complex for branch length optimization in a tree-wide Newton step

### Status

v0: not implemented. v1: not implemented.

---

## Poisson Indel Process (PIP, Bouchard-Cote and Jordan 2012)

A continuous-time Markov process closely related to TKF91 but with a different insertion mechanism. Introduced by <a id="cite-13"></a>[Bouchard-Cote and Jordan 2013](https://doi.org/10.1073/pnas.1220450110) [[13](#ref-13)].

### Mathematical formulation

PIP differs from TKF91 in how insertions are modeled: instead of character births from existing characters (linked birth-death), PIP models insertions as a Poisson process along tree branches. Each inserted character has an independent exponential lifetime (deletion follows the same dynamics as TKF91). This decoupling makes the marginal likelihood computation tractable.

### Key difference from TKF91

- TKF91: insertions are births from existing characters, making positions dependent. Marginal likelihood is $O(L^N)$.
- PIP: insertions are a Poisson process on branches, independent of existing characters. Marginal likelihood is linear in the number of taxa.

The tradeoff: PIP's equilibrium sequence length distribution is Poisson (vs geometric for TKF91), but PIP gains tractability for Bayesian inference.

### Software

- ARPIP: ancestral reconstruction under PIP (<a id="cite-14"></a>[Jowkar et al. 2022](https://doi.org/10.1093/sysbio/syac050) [[14](#ref-14)])
- ProPIP: progressive alignment under PIP (<a id="cite-15"></a>[Maiolo et al. 2018](https://doi.org/10.1186/s12859-018-2357-1) [[15](#ref-15)])

### Practical considerations for branch length optimization

PIP's linear-time marginal likelihood is attractive compared to TKF91's exponential cost, but PIP still requires restructuring the likelihood computation to handle alignment uncertainty. For the branch-length-prevents-zero use case, PIP is over-engineered: a Poisson count model achieves the same practical effect with far less implementation complexity.

### Status

v0: not implemented. v1: not implemented.

---

## General Geometric Indel model (GGI)

The natural generalization of TKF91 to geometrically distributed indel lengths. Described in the review by <a id="cite-16"></a>[Redelings 2024](https://doi.org/10.1093/molbev/msae177) [[16](#ref-16)]. Parameters: insertion rate $\lambda$, deletion rate $\mu$, mean insertion length $\bar{X}$, mean deletion length $\bar{Y}$. TKF91 is the special case $\bar{X} = \bar{Y} = 1$. No exact analytical solution exists; all known approaches use numerical approximation (ODE integration, finite-state machine coarse-graining, pair HMM guesses).

### Status

v0: not implemented. v1: not implemented.

---

## Cumulative Indel Model (De Maio 2020)

Approximates realistic indel dynamics via differential equations with adaptive banding. Introduced by <a id="cite-17"></a>[De Maio 2020](https://doi.org/10.1093/sysbio/syaa050) [[17](#ref-17)]. 18 citations.

### Status

v0: not implemented. v1: not implemented.

---

## Indel-aware parsimony (indelMaP)

A parsimony criterion that treats insertions and deletions as separate events with affine gap penalties for long indels. Introduced by <a id="cite-18"></a>[Iglhaut et al. 2024](https://doi.org/10.1093/molbev/msae109) [[18](#ref-18)]. Separates insertions from deletions on the tree, avoiding the conflation inherent in standard gap-as-missing-data treatment.

### Status

v0: not implemented. v1: not implemented (Fitch parsimony uses majority rule for gap/non-gap at internal nodes).

---

## Binary encoding of indels

Encodes gap presence/absence as a two-state character (0/1) at each alignment position, then analyzes under a two-state substitution model (CFN or Cavender-Farris-Neyman). Mentioned in the "Models of DNA evolution" Wikipedia article.

This is a simplistic approach that ignores indel length, treats each gapped position independently (not grouping contiguous gaps into events), and loses information about insertion vs deletion direction.

### Status

v0: not implemented. v1: not implemented.

---

## Comparison for branch length optimization

| Model                | Parameters           | Indel lengths | Ins/del distinction | Computational cost | Implementation effort | Biological realism |
| :------------------- | :------------------- | :------------ | :------------------ | :----------------- | :-------------------- | :----------------- |
| Gaps as missing data | 0                    | N/A           | No                  | None               | None                  | Ignores indels     |
| Poisson count (v1)   | 1 ($\mu$)            | Equal weight  | No                  | Negligible         | Minimal               | Low (count only)   |
| Affine gap penalty   | 2 ($d$, $e$)         | Geometric     | No                  | Low                | Low                   | Low (heuristic)    |
| Binary encoding      | 0                    | Per-position  | No                  | Low                | Low                   | Very low           |
| TKF91                | 2 ($\lambda$, $\mu$) | Single only   | Yes                 | $O(L^N)$ exact     | High                  | Moderate           |
| TKF92                | 3                    | Geometric     | Yes                 | $O(L^2)$ pairwise  | High                  | Moderate           |
| PIP                  | 2 ($\lambda$, $\mu$) | Single only   | Yes                 | Linear in taxa     | High                  | Moderate           |
| Long indel (SID)     | 4                    | Geometric     | Yes                 | Numerical ODE      | Very high             | High               |
| GGI                  | 4                    | Geometric     | Yes                 | Numerical ODE      | Very high             | High               |

For TreeTime's use case (rapid viral phylodynamics with short alignments and rare indels), the Poisson count model provides the best cost/benefit ratio. The primary goal is preventing zero-length assignment on indel-only branches, not reconstructing the indel process or estimating indel-specific evolutionary parameters. More sophisticated models (TKF91, PIP) would add substantial implementation and computational cost for marginal benefit on typical viral datasets.

For datasets with frequent, long indels (bacterial genomics, gene family evolution), the Poisson count model underweights long indels (a 100-base deletion counts the same as a 1-base deletion). A length-aware Poisson model or an affine-penalty-inspired weighting would be more appropriate. See [proposal](../port-proposals/optimize-indel-model-alternatives.md).

---

## File Index

| File                                                                                                                               | Algorithms                                    |
| :--------------------------------------------------------------------------------------------------------------------------------- | :-------------------------------------------- |
| [`packages/treetime/src/commands/optimize/optimize_indel.rs`](../../packages/treetime/src/commands/optimize/optimize_indel.rs)     | Poisson indel log-likelihood, rate estimation |
| [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs) | Newton integration of indel contribution      |

---

## References

1. <a id="ref-1"></a> Warnow, Tandy. 2012. "Standard Maximum Likelihood Analyses of Alignments with Gaps Can Be Statistically Inconsistent." _PLOS Currents Tree of Life_. https://doi.org/10.1371/currents.rrn1308 [↩](#cite-1)
2. <a id="ref-2"></a> Gotoh, Osamu. 1982. "An Improved Algorithm for Matching Biological Sequences." _Journal of Molecular Biology_ 162(3):705-708. https://doi.org/10.1016/0022-2836(82)90398-9 [↩](#cite-2)
3. <a id="ref-3"></a> Rivas, Elena. 2005. "Evolutionary Models for Insertions and Deletions in a Probabilistic Modeling Framework." _BMC Bioinformatics_ 6:63. https://doi.org/10.1186/1471-2105-6-63 [↩](#cite-3)
4. <a id="ref-4"></a> Qian, Bin, and Richard A. Goldstein. 2001. "Distribution of Indel Lengths." _Proteins: Structure, Function, and Genetics_ 45(1):102-104. https://doi.org/10.1002/prot.1129 [↩](#cite-4)
5. <a id="ref-5"></a> Cartwright, Reed A. 2009. "Problems and Solutions for Estimating Indel Rates and Length Distributions." _Molecular Biology and Evolution_ 26(2):473-480. https://doi.org/10.1093/molbev/msn275 [↩](#cite-5)
6. <a id="ref-6"></a> Wygoda, Elya, Gil Noll, Haim Ashkenazy, Oren Haritan, Jasper Leal, Jotun Hein, Gerton Lunter, and Tal Pupko. 2024. "Statistical Framework to Determine Indel-Length Distribution." _Bioinformatics_ 40(2):btae043. https://doi.org/10.1093/bioinformatics/btae043 [↩](#cite-6)
7. <a id="ref-7"></a> Thorne, Jeffrey L., Hirohisa Kishino, and Joseph Felsenstein. 1991. "An Evolutionary Model for Maximum Likelihood Alignment of DNA Sequences." _Journal of Molecular Evolution_ 33(2):114-124. https://doi.org/10.1007/BF02193625 [↩](#cite-7)
8. <a id="ref-8"></a> Holmes, Ian. 2005. "Using Evolutionary EM to Estimate Indel Rates." _Bioinformatics_ 21(10):2294-2300. https://doi.org/10.1093/bioinformatics/bti177 [↩](#cite-8)
9. <a id="ref-9"></a> Rivas, Elena. 2008. "Probabilistic Phylogenetic Inference with Insertions and Deletions." _PLOS Computational Biology_ 4(9):e1000172. https://doi.org/10.1371/journal.pcbi.1000172 [↩](#cite-9)
10. <a id="ref-10"></a> Thorne, Jeffrey L., Hirohisa Kishino, and Joseph Felsenstein. 1992. "Inching toward Reality: An Improved Likelihood Model of Sequence Evolution." _Journal of Molecular Evolution_ 34(1):3-16. https://doi.org/10.1007/BF00163848 [↩](#cite-10)
11. <a id="ref-11"></a> Miklos, Istvan, Gerton A. Lunter, and Ian Holmes. 2004. "A 'Long Indel' Model for Evolutionary Sequence Alignment." _Molecular Biology and Evolution_ 21(3):529-540. https://doi.org/10.1093/molbev/msh043 [↩](#cite-11)
12. <a id="ref-12"></a> Holmes, Ian. 2020. "A Model of Indel Evolution by Finite-State, Continuous-Time Machines." _Genetics_ 215(4):1187-1204. https://doi.org/10.1534/genetics.120.303630 [↩](#cite-12)
13. <a id="ref-13"></a> Bouchard-Cote, Alexandre, and Michael I. Jordan. 2013. "Evolutionary Inference via the Poisson Indel Process." _Proceedings of the National Academy of Sciences_ 110(4):1160-1166. https://doi.org/10.1073/pnas.1220450110 [↩](#cite-13)
14. <a id="ref-14"></a> Jowkar, Gholamhossein, Jitka Pecerska, Massimo Maiolo, Manuel Gil, and Maria Anisimova. 2022. "ARPIP: Ancestral Sequence Reconstruction with Insertions and Deletions under the Poisson Indel Process." _Systematic Biology_ 72(2):460-472. https://doi.org/10.1093/sysbio/syac050 [↩](#cite-14)
15. <a id="ref-15"></a> Maiolo, Massimo, Xiaolei Zhang, Manuel Gil, and Maria Anisimova. 2018. "Progressive Multiple Sequence Alignment with Indel Evolution." _BMC Bioinformatics_ 19:331. https://doi.org/10.1186/s12859-018-2357-1 [↩](#cite-15)
16. <a id="ref-16"></a> Redelings, Benjamin D. 2024. "Insertions and Deletions: Computational Methods, Evolutionary Dynamics, and Biological Applications." _Molecular Biology and Evolution_ 41(9):msae177. https://doi.org/10.1093/molbev/msae177 [↩](#cite-16)
17. <a id="ref-17"></a> De Maio, Nicola. 2020. "The Cumulative Indel Model: Fast and Accurate Statistical Evolutionary Alignment." _Systematic Biology_ 70(2):236-257. https://doi.org/10.1093/sysbio/syaa050 [↩](#cite-17)
18. <a id="ref-18"></a> Iglhaut, Fynn, Jitka Pecerska, Manuel Gil, and Maria Anisimova. 2024. "Please Mind the Gap: Indel-Aware Parsimony for Fast and Accurate Ancestral Sequence Reconstruction." _Molecular Biology and Evolution_ 41(7):msae109. https://doi.org/10.1093/molbev/msae109 [↩](#cite-18)
