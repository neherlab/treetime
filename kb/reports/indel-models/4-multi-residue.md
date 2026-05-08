# Chapter 4: Multi-residue extensions and approximations

[Back to index](README.md) | Previous: [Chapter 3: Single-residue birth-death models](3-single-residue.md) | Next: [Chapter 5: Empirical models and parsimony](5-empirical.md)

## TKF92 (Thorne-Kishino-Felsenstein 1992)

Extension of TKF91 to multi-residue indels via indivisible "fragments." Introduced by <a id="cite-14"></a>[Thorne, Kishino, and Felsenstein 1992](https://doi.org/10.1007/BF00163848) [[14](#ref-14)].

### Biological model

Sequences are composed of fragments: contiguous blocks of linked residues. An entire fragment is inserted or deleted as a unit. Fragment lengths follow a geometric distribution with parameter $r$ (probability of extending by one residue). TKF91 is the special case $r = 0$ (all fragments have length 1).

The indel process operates on fragments, not individual residues: an insertion creates a new fragment adjacent to an existing one, a deletion removes an entire fragment. Substitutions within fragments follow the standard site-independent model.

### Parameters

$\lambda$ (fragment insertion rate), $\mu$ (fragment deletion rate), $r$ (geometric fragment extension probability). The constraint $\lambda < \mu$ still applies.

### Computational complexity

Pairwise: $O(L_1 \cdot L_2)$ via pair HMM. Multi-sequence: approximate (no exact polynomial solution for general trees).

### Software implementations

**[BAli-Phy](https://github.com/bredelings/BAli-Phy)** implements TKF2 (`src/imodel/imodel.cc:686-732`) by building the TKF1 pair HMM and applying `fragmentize(Q, e)`, which adds geometric self-loops to the M, G1, and G2 states with probability $e = r$. The `lengthp` method is marked `FIXME - this is wrong` and calls `std::abort()`, indicating the equilibrium length distribution for TKF2 is not fully implemented.

**[rust-phylo](https://github.com/acg-team/rust-phylo)** implements TKF92 via `TKF92IndelModel` with fragment length parameter $r$. In TKF92, contiguous gap/non-gap regions are grouped into blocks (vs TKF91 where each column is its own block). The likelihood factorization uses the same subtree product structure as TKF91 but with modified insertion factors: $\lambda\beta(1-r)/r$ replaces $\lambda\beta$.

**[StatAlign](https://github.com/statalign/statalign)** uses TKF92 as its indel model for Bayesian MCMC statistical alignment.

### Relationship to TKF91

The fragment structure addresses TKF91's main biological limitation (single-residue indels) but introduces an approximation: in reality, fragment boundaries are not fixed, and the geometric fragment length distribution is a simplification. The fragment model also does not account for overlapping indels (an insertion within a previously inserted fragment).

---

## RS07 (Redelings-Suchard 2007)

A pair HMM approximation to the General Geometric Indel (GGI) model. Used as the default indel model in BAli-Phy and Historian. Introduced by <a id="cite-15"></a>[Redelings and Suchard 2007](https://doi.org/10.1186/1471-2148-7-40) [[15](#ref-15)].

### Mathematical formulation

RS07 parameterizes a 5-state pair HMM (Start, Match, Gap-in-1, Gap-in-2, End) with two user-facing parameters: `rate` (indel rate relative to substitution rate) and `mean_length` (expected indel length). The geometric extension probability is $\epsilon = (\text{mean\_length} - 1) / \text{mean\_length}$.

The pair HMM is constructed in three steps (as implemented in BAli-Phy `src/builtins/Alignment.cc:166-245`):

1. Compute the scaled rate $\mu = D / (1 - \epsilon)$ where $D = \text{rate} \times t$
2. Compute indel opening probability $\delta = P_{\text{indel}} / (1 + P_{\text{indel}})$ where $P_{\text{indel}} = 1 - e^{-\mu}$
3. Build transition matrix: $Q(S, M) = 1 - 2\delta$, $Q(S, G1) = Q(S, G2) = \delta$
4. Apply `fragmentize(Q, \epsilon)` to add geometric self-loops to M, G1, G2 states
5. Remove the silent Start state by marginalization

### Difference from TKF91

| Aspect           | RS07                                                         | TKF91                              |
| :--------------- | :----------------------------------------------------------- | :--------------------------------- |
| Parameterization | rate + mean_length                                           | $\lambda$ + $\mu$ (birth-death)    |
| Time dependence  | Linear: $D = \text{rate} \times t$, then Poisson gap opening | Exact continuous-time Markov chain |
| Indel lengths    | Geometric via self-loops                                     | Single residue only                |
| Ins/del symmetry | Symmetric ($\delta$ for both)                                | Asymmetric (different G1 row)      |
| Model class      | Phenomenological approximation to GGI                        | Exact birth-death process          |

RS07 is a practical approximation designed for Bayesian MCMC: it captures the essential features of the GGI model (geometric indel lengths, rate proportional to branch length) without the computational cost of solving the GGI differential equations. The approximation is accurate for short branches and deteriorates for long branches where multiple overlapping indels are common.

### Software

**[BAli-Phy](https://github.com/bredelings/BAli-Phy)**: default indel model. **[Historian](https://github.com/ihh/dart)**: uses RS07 for progressive alignment with MCMC refinement. Reference for Historian: <a id="cite-16"></a>[Holmes 2017](https://doi.org/10.1186/s12859-017-1665-1) [[16](#ref-16)].

---

## General Geometric Indel model (GGI)

The natural generalization of TKF91 to geometrically distributed indel lengths. Parameters: insertion rate $\lambda$, deletion rate $\mu$, mean insertion length $\bar{X}$, mean deletion length $\bar{Y}$. TKF91 is the special case $\bar{X} = \bar{Y} = 1$.

### Why GGI has no exact solution

When deletions remove more than one residue, adjacent residue fates become dependent: a single deletion event simultaneously removes multiple residues, correlating their evolutionary histories. The finite-time gap length distribution becomes a convolution of multiple overlapping indel events and is no longer geometric. No finite-state pair HMM can represent GGI exactly.

Time-reversibility holds iff $\lambda(\bar{X} - 1) = \mu(\bar{Y} - 1)$, which constrains the parameter space.

### Approximations to GGI

The intractability of exact GGI has motivated several approximation strategies:

- TKF92: fragments (see above)
- <a id="cite-25"></a>Knudsen and Miyamoto 2003 [[25](#ref-25)], Redelings and Suchard 2005/2007 (RS07): guessed pair HMM forms matching GGI moments (see above)
- De Maio 2021: moment-matching differential equations for best-fit pair HMM. The Cumulative Indel Model approximates GGI dynamics via ODEs with adaptive banding. <a id="cite-21"></a>[De Maio 2021](https://doi.org/10.1093/sysbio/syaa050) [[21](#ref-21)].
- Holmes 2020: refined ODEs via coarse-graining of pair HMM state spaces. <a id="cite-22"></a>[Holmes 2020](https://doi.org/10.1534/genetics.120.303630) [[22](#ref-22)]. The best known approximation to GGI.

### The Redelings 2024 review

<a id="cite-23"></a>[Redelings et al. 2024](https://doi.org/10.1093/molbev/msae177) [[23](#ref-23)] provides a taxonomy of GGI and its approximations. Key recommendation from the review: indel information should only be used when the alignment was inferred with a phylogeny-aware aligner. Ideally MSA and tree should be co-estimated. Indel rates should be inferred from data, not fixed to defaults. Single-character indel models conflate indel rate and indel length into one parameter, making the parameter uninterpretable as a biological rate.

---

## Long indel model (SID, Miklos-Lunter-Holmes 2004)

A continuous-time Markov model allowing indels of arbitrary length with geometric length distributions. Introduced by <a id="cite-24"></a>[Miklos, Lunter, and Holmes 2004](https://doi.org/10.1093/molbev/msh043) [[24](#ref-24)].

### Mathematical formulation

Extends TKF91 by allowing single insertion and deletion events of any length, with geometric length distributions parameterized separately for insertions and deletions. Parameters: insertion rate $\lambda$, deletion rate $\mu$, mean insertion length $\bar{X}$, mean deletion length $\bar{Y}$. Uses an infinite sequence embedding to make rates position-independent and introduces "rate grammar" notation.

The key algorithmic contribution is the "chop zone" decomposition: the alignment is partitioned into zones where the trajectory likelihood can be computed independently, yielding $O(L^2)$ pairwise complexity.

### Relationship to GGI

The SID model and GGI describe the same underlying process. "SID" refers to the Miklos-Lunter-Holmes formulation and algorithmic approach (chop zones, rate grammars). "GGI" is the broader model class. TKF91 is the special case $\bar{X} = \bar{Y} = 1$ of both.

---

## References

14. <a id="ref-14"></a> Thorne, Jeffrey L., Hirohisa Kishino, and Joseph Felsenstein. 1992. "Inching toward Reality." _J Mol Evol_ 34(1):3-16. https://doi.org/10.1007/BF00163848
15. <a id="ref-15"></a> Redelings, Benjamin D., and Marc A. Suchard. 2007. "Incorporating Indel Information into Phylogeny Estimation." _BMC Evol Biol_ 7(1):40. https://doi.org/10.1186/1471-2148-7-40
16. <a id="ref-16"></a> Holmes, Ian. 2017. "Solving the Master Equation for Indels." _BMC Bioinformatics_ 18:255. https://doi.org/10.1186/s12859-017-1665-1
17. <a id="ref-21"></a> De Maio, Nicola. 2021. "The Cumulative Indel Model." _Syst Biol_ 70(2):236-257. https://doi.org/10.1093/sysbio/syaa050
18. <a id="ref-22"></a> Holmes, Ian. 2020. "A Model of Indel Evolution by Finite-State, Continuous-Time Machines." _Genetics_ 215(4):1187-1204. https://doi.org/10.1534/genetics.120.303630
19. <a id="ref-23"></a> Redelings, Benjamin D., Ian Holmes, Gerton Lunter, Tal Pupko, and Maria Anisimova. 2024. "Insertions and Deletions." _MBE_ 41(9):msae177. https://doi.org/10.1093/molbev/msae177
20. <a id="ref-24"></a> Miklos, Istvan, et al. 2004. "A 'Long Indel' Model." _MBE_ 21(3):529-540. https://doi.org/10.1093/molbev/msh043

21. <a id="ref-25"></a> Knudsen, Bjarne, and Michael M. Miyamoto. 2003. "Sequence Alignments and Pair Hidden Markov Models Using Evolutionary History." _J Mol Biol_ 333(2):453-460. https://doi.org/10.1016/j.jmb.2003.08.015

See [consolidated references](references.md) for the complete bibliography.
