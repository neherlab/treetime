# Chapter 3: Single-residue birth-death models

[Back to index](README.md) | Previous: [Chapter 2: Gap treatment and counting models](2-gap-treatment.md) | Next: [Chapter 4: Multi-residue extensions and approximations](4-multi-residue.md)

## TKF91 (Thorne-Kishino-Felsenstein 1991)

The first continuous-time Markov chain model jointly modeling substitutions and indels as a single evolutionary process. Introduced by <a id="cite-7"></a>[Thorne, Kishino, and Felsenstein 1991](https://doi.org/10.1007/BF02193625) [[7](#ref-7)].

### Biological model

Sequence evolution is modeled as a linear birth-death process on a chain of characters, anchored by an <a id="gloss-use-6"></a>"immortal link" <sup>[6](#gloss-6)</sup> at the left end (preventing the entire sequence from being deleted):

- Each existing character can spawn a new adjacent character (insertion) at rate $\lambda$
- Each existing character can be removed (deletion) at rate $\mu$
- Constraint: $\lambda < \mu$ for stationarity (finite expected sequence length)
- Substitutions at each site follow a standard model (JC69, HKY85, GTR, etc.) independently of the indel process

The equilibrium sequence length distribution is geometric with parameter $\sigma = \lambda/\mu$: $P(L = l) = (1 - \sigma)\sigma^l$.

### Transition probabilities

The finite-time transition probabilities for a branch of length $t$ are analytically tractable. The key intermediate quantity:

$$\beta(t) = \frac{1 - e^{(\lambda - \mu)t}}{\mu - \lambda e^{(\lambda - \mu)t}}$$

From $\beta$, the event probabilities are:

$$p_1(t) = e^{-\mu t}(1 - \lambda \beta) \quad \text{(homolog survival: ancestor residue present in descendant)}$$
$$p_0'(t) = \mu \beta \quad \text{(deletion: ancestor residue absent in descendant)}$$
$$p_1''(t) = 1 - \lambda \beta \quad \text{(no insertion to right of immortal link)}$$

The insertion probability (a new residue appears in the descendant to the right of a surviving ancestor residue) is $\lambda\beta$. The deletion-then-insertion correction factor $\eta$ accounts for the case where an ancestor residue is deleted but a new residue is inserted in its place. The intermediate product $n_1$ is:

$$n_1 = (1 - e^{-\mu t} - \mu\beta)(1 - \lambda\beta)$$
$$\eta = \ln(n_1) - \ln(\lambda\beta) - \ln(\mu\beta)$$

### Pair HMM representation

TKF91 can be represented as a three-state <a id="gloss-use-3"></a>pair HMM <sup>[3](#gloss-3)</sup> with states Match (M), Insert (G1), and Delete (G2). The transition matrix, as implemented in BAli-Phy (`src/imodel/imodel.cc:579-638`):

| From \ To  | M                                               | G1             | G2                                                  | End                                       |
| :--------- | :---------------------------------------------- | :------------- | :-------------------------------------------------- | :---------------------------------------- |
| Start/M/G2 | $(1-\lambda\beta)\frac{\lambda}{\mu}e^{-\mu t}$ | $\lambda\beta$ | $(1-\lambda\beta)\frac{\lambda}{\mu}(1-e^{-\mu t})$ | $(1-\lambda\beta)(1-\frac{\lambda}{\mu})$ |
| G1         | $\frac{\lambda\beta e^{-\mu t}}{1-e^{-\mu t}}$  | $\lambda\beta$ | $\frac{1-e^{-\mu t}-\mu\beta}{1-e^{-\mu t}}$        | $\frac{(\mu-\lambda)\beta}{1-e^{-\mu t}}$ |

The Start, Match, and G2 rows are identical (the probability of the next column type depends only on whether the ancestor residue survived, not on how it was created). The G1 row is different because it conditions on the ancestor residue having been inserted (not a continuation of a surviving lineage).

The pairwise alignment likelihood under TKF91 is the forward probability of this pair HMM, computed in $O(L_1 \cdot L_2)$ time.

### Multi-sequence complexity

For $N$ sequences on a tree, the exact marginal likelihood (summing over all possible alignments) requires an $N$-dimensional DP hypercube: $O(L^N)$, which is exponential in the number of taxa. This is equivalent to Sankoff's <a id="cite-8"></a>[1975](https://doi.org/10.1137/0128004) [[8](#ref-8)] simultaneous alignment and phylogeny problem, which is NP-complete for the optimal-scoring path.

With a fixed alignment (treating the alignment as observed data, not marginalizing over it), the likelihood computation reduces to an $O(L \cdot N)$ postorder traversal. <a id="cite-9"></a>[Rivas 2008](https://doi.org/10.1371/journal.pcbi.1000172) [[9](#ref-9)] showed how to extend <a id="gloss-use-5"></a>Felsenstein's peeling <sup>[5](#gloss-5)</sup> algorithm to handle single-residue indel events per alignment column under a non-reversible birth-death model.

<a id="cite-10"></a>[Westesson et al. 2012](https://doi.org/10.1371/journal.pone.0034572) [[10](#ref-10)] showed the equivalence of TKF91 to a pair HMM and generalized Felsenstein's peeling from single sites to entire sequences via phylogenetic <a id="gloss-use-4"></a>transducers <sup>[4](#gloss-4)</sup>, reducing the multi-sequence complexity to $O(L^2 N)$ under MCMC.

### Software implementations

**[BAli-Phy](https://github.com/bredelings/BAli-Phy)** (C++). Bayesian joint alignment and phylogeny inference via MCMC. Implements TKF91 as `TKF1` in `src/imodel/imodel.cc:579-638`. Parameters: $\lambda = \exp(\text{parameter}[0])$ (log-scale), mean sequence length from which $\mu = \lambda/\sigma$ is derived. Builds the full pair HMM transition matrix. Branch lengths are jointly estimated with alignment and indel parameters via MCMC, not Newton-Raphson. Reference: <a id="cite-11"></a>[Redelings and Suchard 2005](https://doi.org/10.1080/10635150590947041) [[11](#ref-11)].

**[BEAST](https://github.com/beast-dev/beast-mcmc)** (Java). `dr.oldevomodel.indel.TKF91Likelihood` wraps `HomologyRecursion` (680 lines). Uses `BFloat` arbitrary-precision arithmetic to avoid underflow in the DP table. Parameters: `lambda = deathRate * lengthDistributionValue`, `mu = deathRate`. Precomputes per-node arrays ($\beta$, $\mu\beta$, $e^{-\mu t}(1-\lambda\beta)$, etc.) for the tree recursion.

**[rust-phylo](https://github.com/acg-team/rust-phylo)** (Rust). `TKF91IndelModel` in `phylo/src/tkf_model/tkf91.rs`. Uses `nalgebra` `DMatrix` for intermediate storage. The TKF likelihood factorizes over alignment blocks via subtree products: root term $\ln(1-\lambda/\mu)$, plus immortal link terms $\sum \ln(1-\lambda\beta)$, plus per-block event factors accumulated in a postorder traversal. Branch length optimization uses Brent's method (`argmin::BrentOpt`) with the combined TKF+substitution cost function. Supports TKF92 extension via a fragment length parameter $r$.

**[CRAN/TKF](https://github.com/cran/TKF)** (C/R). Pairwise distance estimation. Three DP matrices (gap-in-B, match, gap-in-A) in log space. Uses GSL Brent minimizer for 1D or Nelder-Mead/BFGS for 2D ML estimation of $(\lambda, \mu)$ or $(t, \mu)$.

**[StatAlign](https://github.com/statalign/statalign)** (Java). Bayesian MCMC statistical alignment. Uses TKF92 (fragment model). Reference: <a id="cite-12"></a>[Novak et al. 2008](https://doi.org/10.1093/bioinformatics/btn457) [[12](#ref-12)].

### Applicability to TreeTime branch length optimization

TKF91 distinguishes insertions from deletions (separate $\lambda$ and $\mu$), which is more biologically realistic than the Poisson count model's single rate. The pair HMM formalism integrates naturally with tree likelihood computation. The main barriers to adoption in TreeTime:

1. **Architectural mismatch.** TreeTime's per-edge optimizer uses eigendecomposition-based coefficients cached once per edge. TKF91 requires pair HMM DP per edge, which has different data flow (sequence-level, not site-level).
2. **Parameter estimation.** TKF91 adds two parameters ($\lambda$, $\mu$) that must be estimated jointly with branch lengths. The current Poisson model estimates one parameter ($\mu_{\text{indel}}$) analytically.
3. **Single-residue limitation.** TKF91 allows only single-residue indels. This is less realistic than the Poisson count model for datasets with multi-residue indels (the Poisson model at least counts each multi-residue event as one event; TKF91 cannot represent them at all without the TKF92 extension).

### Additional references

- <a id="cite-7b"></a>[Thorne, Kishino, and Felsenstein 1991](https://doi.org/10.1007/BF02193625) [[7](#ref-7)] (TKF91 original paper)
- <a id="cite-13"></a>[Holmes 2005](https://doi.org/10.1093/bioinformatics/bti177) [[13](#ref-13)] (EM estimation of TKF91 insertion/deletion rates)

---

## Poisson Indel Process (PIP, Bouchard-Cote and Jordan 2013)

A continuous-time Markov process that achieves tractable marginal likelihood computation by decoupling insertions from existing characters. Introduced by <a id="cite-17"></a>[Bouchard-Cote and Jordan 2013](https://doi.org/10.1073/pnas.1220450110) [[17](#ref-17)].

### Biological model

PIP differs from TKF91 in how insertions are modeled. In TKF91, each existing character can spawn a new adjacent character (linked birth process), making the fates of adjacent positions dependent. In PIP, insertions are a Poisson process along tree branches: new characters appear at rate $\lambda$ per unit branch length, independent of existing characters. Each character (whether original or inserted) has an independent exponential lifetime with rate $\mu$ (deletion follows the same dynamics as TKF91).

The decoupling of insertions from existing characters is the key theoretical contribution: it makes the marginal likelihood (summing over all possible alignments) computable in time linear in the number of taxa, vs $O(L^N)$ for TKF91.

### Tradeoff vs TKF91

The price of tractability is that PIP's equilibrium sequence length distribution is Poisson (vs geometric for TKF91). The Poisson distribution has lighter tails than geometric, which means PIP predicts less length variation at equilibrium. In practice, this matters less than the computational advantage, especially for Bayesian inference where marginalizing over alignments is essential.

### Software

**[ProPIP](https://github.com/acg-team/ProPIP)**: ML progressive alignment under PIP. The first polynomial-time progressive aligner with a rigorous indel model. <a id="cite-18"></a>[Maiolo et al. 2018](https://doi.org/10.1186/s12859-018-2357-1) [[18](#ref-18)]. Estimates indel rates from data. Supports Gamma rate heterogeneity.

**[ARPIP](https://github.com/acg-team/bpp-arpip)**: ancestral sequence reconstruction under PIP. <a id="cite-19"></a>[Jowkar et al. 2022](https://doi.org/10.1093/sysbio/syac050) [[19](#ref-19)]. Two-step approach: find most probable indel points, then reconstruct on the pruned subtree. <a id="cite-20"></a>[Jowkar et al. 2024](https://doi.org/10.1186/s12859-024-05986-1) [[20](#ref-20)] showed that the single-character indel assumption still captures long indel patterns on mammalian protein orthologs.

**[rust-phylo](https://github.com/acg-team/rust-phylo)** implements PIP alongside TKF91/TKF92.

### Applicability to TreeTime

PIP's linear-time marginal likelihood is attractive compared to TKF91's exponential cost. But PIP still requires alignment-aware likelihood computation (tracking which characters are present on which branches), not just a per-edge scalar. The integration would change TreeTime's optimization architecture: the current per-edge eigendecomposition coefficients cache would need to be replaced or augmented with PIP's column-based likelihood. For the branch-length-prevents-zero use case, the Poisson count model achieves the same practical effect with negligible implementation cost.

---

## Glossary

1. <a id="gloss-6"></a> **Immortal link.** A TKF91 concept: a permanent anchor at the left end of the sequence that cannot be deleted. [↩](#gloss-use-6)
2. <a id="gloss-3"></a> **Pair HMM.** A hidden Markov model generating two sequences simultaneously, with Match/Insert/Delete states. [↩](#gloss-use-3)
3. <a id="gloss-5"></a> **Felsenstein peeling.** Postorder tree traversal for phylogenetic likelihoods ([Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[29](#ref-29)]). [↩](#gloss-use-5)
4. <a id="gloss-4"></a> **Phylogenetic transducer.** Finite-state machine transforming ancestor into descendant sequence ([Bradley and Holmes 2007](https://doi.org/10.1093/bioinformatics/btm402) [[28](#ref-28)]). [↩](#gloss-use-4)

## References

7. <a id="ref-7"></a> Thorne, Jeffrey L., Hirohisa Kishino, and Joseph Felsenstein. 1991. "An Evolutionary Model for Maximum Likelihood Alignment of DNA Sequences." _J Mol Evol_ 33(2):114-124. https://doi.org/10.1007/BF02193625
8. <a id="ref-8"></a> Sankoff, David. 1975. "Minimal Mutation Trees of Sequences." _SIAM J Appl Math_ 28(1):35-42. https://doi.org/10.1137/0128004
9. <a id="ref-9"></a> Rivas, Elena. 2008. "Probabilistic Phylogenetic Inference with Insertions and Deletions." _PLOS Comp Biol_ 4(9):e1000172. https://doi.org/10.1371/journal.pcbi.1000172
10. <a id="ref-10"></a> Westesson, Oscar, et al. 2012. "Accurate Reconstruction of Insertion-Deletion Histories." _PLOS ONE_ 7(4):e34572. https://doi.org/10.1371/journal.pone.0034572
11. <a id="ref-11"></a> Redelings, Benjamin D., and Marc A. Suchard. 2005. "Joint Bayesian Estimation of Alignment and Phylogeny." _Syst Biol_ 54(3):401-418. https://doi.org/10.1080/10635150590947041
12. <a id="ref-12"></a> Novak, Adam, et al. 2008. "StatAlign." _Bioinformatics_ 24(20):2403-2404. https://doi.org/10.1093/bioinformatics/btn457
13. <a id="ref-13"></a> Holmes, Ian. 2005. "Using Evolutionary EM to Estimate Indel Rates." _Bioinformatics_ 21(10):2294-2300. https://doi.org/10.1093/bioinformatics/bti177
14. <a id="ref-17"></a> Bouchard-Cote, Alexandre, and Michael I. Jordan. 2013. "Evolutionary Inference via the Poisson Indel Process." _PNAS_ 110(4):1160-1166. https://doi.org/10.1073/pnas.1220450110
15. <a id="ref-18"></a> Maiolo, Massimo, et al. 2018. "Progressive Multiple Sequence Alignment with Indel Evolution." _BMC Bioinformatics_ 19:331. https://doi.org/10.1186/s12859-018-2357-1
16. <a id="ref-19"></a> Jowkar, Gholamhossein, et al. 2022. "ARPIP." _Syst Biol_ 72(2):460-472. https://doi.org/10.1093/sysbio/syac050
17. <a id="ref-20"></a> Jowkar, Gholamhossein, et al. 2024. "Single-Character Insertion-Deletion Model Preserves Long Indels." _BMC Bioinformatics_ 25:273. https://doi.org/10.1186/s12859-024-05986-1
18. <a id="ref-28"></a> Bradley, Robert K., and Ian Holmes. 2007. "Transducers." _Bioinformatics_ 23(23):3258-3262. https://doi.org/10.1093/bioinformatics/btm402
19. <a id="ref-29"></a> Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences." _J Mol Evol_ 17(6):368-376. https://doi.org/10.1007/BF01734359

See [consolidated references](references.md) for the complete bibliography.
