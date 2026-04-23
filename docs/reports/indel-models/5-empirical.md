# Chapter 5: Empirical models and parsimony

[Back to index](_index.md) | Previous: [Chapter 4: Multi-residue extensions](4-multi-residue.md) | Next: [Chapter 6: Comparative analysis](6-comparison.md)

## SIM/RIM (Loewenthal et al. 2021)

An empirical indel rate model that separates insertion from deletion dynamics and uses power-law (<a id="gloss-use-1"></a>Zipfian <sup>[†1](#gloss-1)</sup>) length distributions. Introduced by <a id="cite-25"></a>[Loewenthal et al. 2021](https://doi.org/10.1093/molbev/msab266) [[25](#ref-25)].

### Models

**SIM (Simple Indel Model)** - 3 parameters:

- $R_{ID}$: indel-to-substitution-rate ratio (insertion rate = deletion rate)
- $A_{ID}$: Zipfian length distribution exponent (same for insertions and deletions)
- $RL$: root sequence length

**RIM (Rich Indel Model)** - 5 parameters:

- $R_I$: insertion-to-substitution-rate ratio
- $R_D$: deletion-to-substitution-rate ratio
- $A_I$: Zipfian length distribution exponent for insertions
- $A_D$: Zipfian length distribution exponent for deletions
- $RL$: root sequence length

Both models parameterize indel rates relative to the substitution rate, not as absolute rates. The truncated Zipfian distribution assigns probability $P(k) \propto k^{-a}$ to indel length $k$, with maximum indel size of 50 amino acids.

### Key empirical findings

Analysis of 4,823 protein datasets across 15 taxonomic groups:

- 35% of datasets favor RIM over SIM (separate ins/del parameters)
- Among RIM datasets: **74% have deletion rate > insertion rate** ($R_D > R_I$)
- The rate asymmetry ($R_D/R_I \approx 2$ in Drosophilidae) is the dominant signal; length asymmetry is weak (only 55% of RIM datasets show $A_D > A_I$)
- Exception: Saccharomycetaceae coding genes show $R_I > R_D$ in 56% of RIM datasets (yeast introns reverse this: 89.5% show $R_D > R_I$)

### Parameter estimation

Inference uses <a id="gloss-use-2"></a>ABC <sup>[†2](#gloss-2)</sup> (Approximate Bayesian Computation), not maximum likelihood. The pipeline simulates 100,000 MSAs by drawing parameters from priors and evolving indels along the input phylogeny via the <a id="gloss-use-3"></a>Gillespie algorithm <sup>[†3](#gloss-3)</sup>, computes 27 summary statistics for each, accepts the closest 0.1%, and averages accepted parameters. LASSO regression corrects for alignment-induced bias in summary statistics. Running time: ~328 minutes per dataset.

### Implications for TreeTime

The SIM/RIM findings validate the direction of TreeTime v1's Poisson indel model but highlight its limitations:

1. **Separate rates matter.** The empirical $R_D/R_I \approx 2$ ratio means the symmetric Poisson model's single rate $\mu_{\text{indel}}$ is a compromise between two distinct processes.
2. **Zipfian lengths matter.** The equal-weight assumption (all events count as 1) discards the signal in indel lengths. A length-weighted Poisson model would better capture the contribution of long indels.
3. **Rate-relative parameterization.** SIM/RIM parameterize indel rates relative to the substitution rate. This is more interpretable than the absolute rate $\mu_{\text{indel}}$ used in v1, because it scales naturally with evolutionary distance.

These observations motivate the [extensions E1 and E2 in the alternatives proposal](../../port-proposals/optimize-indel-model-alternatives.md).

### Software

**[SpartaABC](https://github.com/gilloe/SpartaABC)** (C++ + Python). Webserver: https://spartaabc.tau.ac.il/. Uses the Gillespie algorithm for simulation and scikit-learn for LASSO regression.

---

## Indel-aware parsimony (indelMaP, Iglhaut et al. 2024)

A parsimony criterion that treats insertions and deletions as separate evolutionary events with affine gap costs for long indels. Introduced by <a id="cite-26"></a>[Iglhaut et al. 2024](https://doi.org/10.1093/molbev/msae109) [[26](#ref-26)].

### Approach

Uses the <a id="gloss-use-4"></a>Dollo principle <sup>[†4](#gloss-4)</sup>: a deleted character cannot reappear (no re-insertion). Insertions and deletions are inferred separately on the tree using affine gap penalties (opening + extension costs). This avoids the conflation inherent in standard gap-as-missing-data treatment, where it is impossible to distinguish an ancestral insertion from an ancestral deletion.

### Practical performance

According to the <a id="cite-23"></a>[Redelings et al. 2024](https://doi.org/10.1093/molbev/msae177) [[23](#ref-23)] review, indelMaP outperforms competitors on large densely sampled datasets, where the parsimony assumption (few changes per branch) holds well. For sparse datasets with long branches, probabilistic models (TKF91, RS07, PIP) are more accurate.

### Software

**[indelMaP](https://github.com/acg-team/indelMaP)** (Python). Also available in Rust: [rust-indelMaP](https://github.com/acg-team/rust-indelMaP).

---

## Binary encoding of indels

Encodes gap presence/absence as a two-state character (0/1) at each alignment position, then analyzes under a two-state substitution model (Cavender-Farris-Neyman). Used by [FastML](http://fastml.tau.ac.il/) for indel reconstruction (web server only, no public repository).

This approach ignores indel length (treating each gapped position independently, not grouping contiguous gaps into events), loses information about insertion vs deletion direction, and cannot model overlapping indels. It is the simplest approach that treats indels as evolutionary signal rather than missing data, but the per-position independence assumption is biologically wrong (indel events affect contiguous blocks).

---

## Glossary

1. <a id="gloss-1"></a> **Zipf distribution (power law).** A discrete probability distribution where the probability of observing value $k$ is $P(k) \propto k^{-a}$ for exponent $a > 1$. Larger $a$ means shorter indels dominate; smaller $a$ means long indels are more frequent. Empirical indel length data fits Zipf better than geometric in most datasets ([Cartwright 2009](https://doi.org/10.1093/molbev/msn275) [[5](#ref-5)]; [Wygoda et al. 2024](https://doi.org/10.1093/bioinformatics/btae043) [[6](#ref-6)]). [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **ABC (Approximate Bayesian Computation).** An inference method that bypasses likelihood computation by simulating data from prior parameter draws, comparing simulated and observed summary statistics, and accepting parameter values that produce close matches. [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **Gillespie algorithm.** A stochastic simulation algorithm for continuous-time Markov chains (<a id="cite-30"></a>[Gillespie 1977](https://doi.org/10.1021/j100540a008) [[30](#ref-30)]). Generates exact sample paths by drawing exponentially distributed waiting times between events and choosing the next event type proportional to its rate. [↩](#gloss-use-3)
4. <a id="gloss-4"></a> **Dollo principle (Dollo parsimony).** The assumption that a complex feature, once lost, cannot be regained. Applied to indels: a deleted character cannot be re-inserted at the same position. [↩](#gloss-use-4)

## References

5. <a id="ref-5"></a> Cartwright, Reed A. 2009. "Problems and Solutions for Estimating Indel Rates and Length Distributions." _Molecular Biology and Evolution_ 26(2):473-480. https://doi.org/10.1093/molbev/msn275
6. <a id="ref-6"></a> Wygoda, Elya, et al. 2024. "Statistical Framework to Determine Indel-Length Distribution." _Bioinformatics_ 40(2):btae043. https://doi.org/10.1093/bioinformatics/btae043
7. <a id="ref-23"></a> Redelings, Benjamin D., Ian Holmes, Gerton Lunter, Tal Pupko, and Maria Anisimova. 2024. "Insertions and Deletions: Computational Methods, Evolutionary Dynamics, and Biological Applications." _Molecular Biology and Evolution_ 41(9):msae177. https://doi.org/10.1093/molbev/msae177
8. <a id="ref-25"></a> Loewenthal, Gil, et al. 2021. "A Probabilistic Model for Indel Evolution: Differentiating Insertions from Deletions." _Molecular Biology and Evolution_ 38(12):5769-5781. https://doi.org/10.1093/molbev/msab266
9. <a id="ref-26"></a> Iglhaut, Clara, et al. 2024. "Please Mind the Gap: Indel-Aware Parsimony for Fast and Accurate Ancestral Sequence Reconstruction." _Molecular Biology and Evolution_ 41(7):msae109. https://doi.org/10.1093/molbev/msae109
10. <a id="ref-30"></a> Gillespie, Daniel T. 1977. "Exact Stochastic Simulation of Coupled Chemical Reactions." _Journal of Physical Chemistry_ 81(25):2340-2361. https://doi.org/10.1021/j100540a008

See [consolidated references](references.md) for the complete bibliography across all chapters.
