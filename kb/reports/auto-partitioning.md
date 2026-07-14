# Site-Specific Models and Automatic Partitioning for TreeTime v1

## 1. Problem and motivation

A single <a id="gloss-use-1"></a>GTR <sup>[1](#gloss-1)</sup> substitution model applied uniformly across an alignment ignores the fact that different positions evolve under different selective constraints. This is not just a rate problem. <a id="cite-1a"></a>[Halpern and Bruno 1998](https://doi.org/10.1093/oxfordjournals.molbev.a025995) [[1](#ref-1)] showed that site-to-site differences in amino acid frequencies are the root cause: rate heterogeneity is a consequence of preference heterogeneity, not an independent phenomenon. A model capturing per-site preferences automatically captures much of the rate variation.

When site-specific preferences are ignored, branch lengths are systematically underestimated. The probability of observing the same state at both ends of a branch by chance is $\sum_i (\pi_i^a)^2$, which increases sharply at constrained sites. A model using uniform frequencies misinterprets high identity as close relatedness (<a id="cite-2a"></a>[Puller, Sagulenko, and Neher 2020](https://doi.org/10.1093/ve/veaa066) [[2](#ref-2)], Figs 3-4). <a id="gloss-use-2"></a>Gamma rate variation <sup>[2](#gloss-2)</sup> alone does not fix this because it models rate, not preference (<a id="cite-3"></a>[Yang 1994](https://doi.org/10.1007/BF00178256) [[3](#ref-3)]).

In molecular clock analyses, branch length underestimation creates an apparent time-dependent rate: short branches near the tips appear to evolve faster than long branches near the root, biasing the clock rate downward (<a id="cite-2b"></a>[Puller, Sagulenko, and Neher 2020](https://doi.org/10.1093/ve/veaa066) [[2](#ref-2)]; <a id="cite-4"></a>[Hilton and Bloom 2018](https://doi.org/10.1093/ve/vey033) [[4](#ref-4)]). For HIV, branches longer than ~100 years become inaccurate, and signal disappears beyond ~300 years.

TreeTime v1 currently uses a single GTR model per sequence partition. The goal is to support both site-specific models and automatic partitioning with families of linked GTR models.

---

## 2. Theoretical foundations

### 2.1 Mutation-selection balance

<a id="cite-5"></a>[Bruno 1996](https://doi.org/10.1093/oxfordjournals.molbev.a025583) [[5](#ref-5)] first showed that phylogenetic correlations distort column frequencies in alignments and proposed ML correction. <a id="cite-1b"></a>[Halpern and Bruno 1998](https://doi.org/10.1093/oxfordjournals.molbev.a025995) [[1](#ref-1)] extended this to a full model where position-specific amino acid frequencies are free parameters determined by <a id="gloss-use-3"></a>mutation-selection balance <sup>[3](#gloss-3)</sup>. Their key insight: site-specific frequencies determine the maximum dissimilarity expected for diverged but conserved sequences, and are required for estimating long evolutionary distances.

### 2.2 The site-specific GTR parametrization

The rate matrix at site $a$ is:

$$Q_{ij}^a = \mu^a \, \pi_i^a \, W_{ij} \quad \text{for } i \neq j$$

where $W_{ij}$ is a shared symmetric <a id="gloss-use-4"></a>exchangeability matrix <sup>[4](#gloss-4)</sup>, $\pi_i^a$ is the <a id="gloss-use-5"></a>equilibrium frequency <sup>[5](#gloss-5)</sup> of state $i$ at site $a$, and $\mu^a$ is the overall substitution rate at site $a$. Normalization: $\sum_i \pi_i^a = 1$ and $L^{-1} \sum_a \sum_{i \neq j} W_{ij} \pi_i^a \pi_j^a = 1$.

Sharing $W_{ij}$ across all sites is the key regularizer. Full per-site $W_{ij}^a$ would be unidentifiable for pathogen datasets and computationally prohibitive.

### 2.3 Intuitive understanding

The site-specific model has three kinds of parameters, each with a different physical meaning and a different answer to the question "should this be shared across alignment columns or estimated separately?"

$W_{ij}$ (exchangeability) encodes how easily one nucleotide converts into another. A-to-G transitions are common because purines share a similar molecular structure; A-to-C transversions are rare because the molecules differ more. This reflects nucleotide chemistry, which does not change from one alignment column to the next. Sharing $W$ across all sites is physically motivated and is the key regularizer that keeps the model identifiable with limited data (<a id="cite-2e"></a>[Puller, Sagulenko, and Neher 2020](https://doi.org/10.1093/ve/veaa066) [[2](#ref-2)]).

$\pi_i^a$ (equilibrium frequencies) describes the preferred distribution of nucleotides at a specific position. A structurally constrained position in a protein might strongly prefer adenine, while a position in a hypervariable loop might accept any nucleotide roughly equally. These preferences arise from site-specific fitness landscapes shaped by natural selection (<a id="cite-1c"></a>[Halpern and Bruno 1998](https://doi.org/10.1093/oxfordjournals.molbev.a025995) [[1](#ref-1)]). Different positions have different constraints, so $\pi$ must vary per site or per partition.

$\mu^a$ (overall rate) is a scalar multiplier controlling how fast a position evolves. Third codon positions, where most substitutions are synonymous, evolve faster than first and second positions, where substitutions often change the amino acid. Rate variation across sites was the first form of heterogeneity modeled in phylogenetics (<a id="cite-3c"></a>[Yang 1994](https://doi.org/10.1007/BF00178256) [[3](#ref-3)]).

"Sharing" a parameter across partitions means pooling data from all those partitions to estimate it. Pooling increases the effective sample size and reduces estimation variance, but produces a wrong average if the parameter genuinely differs between partitions. $W$ is safe to pool (chemistry is universal). $\pi$ is risky to pool (different positions have different preferences). $\mu$ should not be pooled (positions evolve at different speeds by definition). This is why the default configuration is: share $W$ broadly, share $\pi$ within groups of functionally similar sites, estimate $\mu$ separately per partition. The [grouping specification](#gtr-family-inference) makes this configurable.

One more sharing axis: branch lengths. Each edge in the tree represents evolutionary distance. With scaled (edge-proportional) branch lengths, all partitions share the same relative branch lengths but each gets a rate multiplier $\mu$. With edge-unlinked branch lengths, each partition gets fully independent branch lengths. Scaled is the standard default (RAxML-NG, IQ-TREE) because partitions on the same tree share the same phylogenetic history; only the pace differs.

### 2.4 Hierarchy of heterogeneity approaches

The following models form a spectrum of increasing complexity:

- Gamma (+G): $k$ rate categories from discretized Gamma($\alpha$), same frequencies for all sites (<a id="cite-3b"></a>[Yang 1994](https://doi.org/10.1007/BF00178256) [[3](#ref-3)])
- <a id="gloss-use-6"></a>FreeRate <sup>[6](#gloss-6)</sup> (+R): $k$ rate categories with freely estimated rates and weights, same frequencies
- Empirical profile mixtures (C10-C60, LG+C60): $k$ fixed frequency profiles from training data, class weights estimated
- CAT (<a id="cite-6"></a>[Lartillot and Philippe 2004](https://doi.org/10.1093/molbev/msh112) [[6](#ref-6)]): infinite mixture via <a id="gloss-use-7"></a>Dirichlet process <sup>[7](#gloss-7)</sup> prior, number of classes adapts to data. Bayesian only. Suppresses long-branch attraction artifacts (<a id="cite-7"></a>[Lartillot, Brinkmann, and Philippe 2007](https://doi.org/10.1186/1471-2148-7-S1-S4) [[7](#ref-7)])
- <a id="gloss-use-8"></a>PMSF <sup>[8](#gloss-8)</sup> (<a id="cite-8a"></a>[Wang et al. 2018](https://doi.org/10.1093/sysbio/syx068) [[8](#ref-8)]): per-site profiles from mixture model posterior means. ML-compatible, fixed during tree search
- Full site-specific (<a id="cite-2c"></a>[Puller, Sagulenko, and Neher 2020](https://doi.org/10.1093/ve/veaa066) [[2](#ref-2)]): per-site $\pi_i^a$ and $\mu^a$, shared $W_{ij}$, estimated from data via EM

Rate-only models are strict subsets: setting $\pi_i^a = \pi_i$ for all $a$ recovers Gamma/FreeRate. Mixture models approximate the site-specific model by clustering sites into groups. PMSF bridges the two: it starts from a finite mixture but produces per-site profiles.

---

## 3. Inference algorithm and convergence

### 3.1 Iterative update rules

<a id="cite-2d"></a>[Puller, Sagulenko, and Neher 2020](https://doi.org/10.1093/ve/veaa066) [[2](#ref-2)] derive multiplicative update rules (Equation 2):

$$\mu^a \leftarrow \mu^a \frac{c + \sum_{ij} n_{ij}^a}{\sum_{ij} \mu^a \pi_i W_{ij} \tau_j^a}$$

$$\pi_i^a \leftarrow \pi_i^a \frac{c + \delta_{is^a} + \sum_{j \neq i} n_{ij}^a}{\pi_i^a (qc + 1) + \mu^a \pi_i^a \sum_j W_{ij} \tau_j^a}$$

$$W_{ij} \leftarrow W_{ij} \frac{\sum_a (n_{ij}^a + n_{ji}^a)}{\sum_a \mu^a W_{ij} (\pi_i^a \tau_j^a + \pi_j^a \tau_i^a)}$$

where $\tau_j^a$ is the time site $a$ spends in state $j$ across the tree, $n_{ij}^a$ is the number of $j \to i$ transitions at site $a$, $s^a$ is the root state, $q$ is the alphabet size, and $c$ is a pseudocount (Dirichlet prior).

### 3.2 Relationship to non-negative matrix factorization

These updates are instances of <a id="gloss-use-9"></a>NMF <sup>[9](#gloss-9)</sup> algorithms (referencing <a id="cite-9"></a>[Lee and Seung 2001](https://proceedings.neurips.cc/paper/2000/hash/f9d1152547c0bde01830b7e8bd60024c-Abstract.html) [[9](#ref-9)]). Each rule is the ratio of observed to expected counts, multiplied by the current estimate. This guarantees positivity, monotone likelihood increase, and convergence to a local maximum.

### 3.3 Accuracy scaling and practical thresholds

Estimation accuracy scales as $\chi^2 \sim (\text{tree length})^{-1}$, limited by Poisson statistics of observable mutations. Practical thresholds:

- Tree length $> 10$: site-specific preferences estimated accurately (most sites mutate multiple times)
- Root-to-tip distance $< 0.3$: naive ancestral reconstruction with a simple model suffices
- Root-to-tip $0.3$ to $1.0$: iterative refinement (model inference + ancestral reconstruction alternation) needed
- Root-to-tip $> 1.0$: tree reconstruction itself becomes unreliable

### 3.4 Iterative refinement loop

1. Reconstruct ancestral states using current model
2. Accumulate per-site $n_{ij}^a$ and $\tau_j^a$ from the reconstruction
3. Update $W_{ij}$, $\pi_i^a$, $\mu^a$ using update rules
4. Re-optimize branch lengths using updated model
5. Repeat until likelihood converges

---

## 4. Identifiability and fundamental limits

### 4.1 The Darwinian uncertainty principle

<a id="cite-10"></a>[Gascuel and Steel 2020](https://doi.org/10.1093/sysbio/syz054) [[10](#ref-10)] prove that root state and rates of state change can each be estimated with high accuracy individually, but not simultaneously from tip data. The optimal global rates for the two tasks have opposite trends.

For TreeTime's use case (tip-dated pathogen sequences), external time calibration partially mitigates this: the molecular clock provides an independent constraint on rates that breaks the rate-pattern trade-off. The risk concentrates on deep branches connecting subtypes or viral lineages.

### 4.2 Practical identifiability problems

From Puller et al. 2020 simulations:

- Skewed frequencies, low rates, and short branches have similar effects on the likelihood
- Using true rates + inferred preferences gives better branch lengths than inferring both jointly
- Even with accurate preference estimates, branch lengths remain underestimated
- Rate misspecification has immediate impact; preference deviations of ~30% are needed before substantial branch length errors appear

---

## 5. Software implementations

### 5.1 Design comparison table

| Axis                      | IQ-TREE                                        | RAxML-NG                                 | PhyloBayes                              | BEAST2                      |
| ------------------------- | ---------------------------------------------- | ---------------------------------------- | --------------------------------------- | --------------------------- |
| $W$ sharing               | Shared across mixture classes                  | Per-partition                            | Shared across all components            | Per-partition               |
| $\pi$ sharing             | Per-class (C10-C60) or per-site (PMSF)         | Per-partition                            | Per-site via Dirichlet process          | Per-partition               |
| $\mu$ sharing             | Per-partition rate scalar                      | Per-partition rate scalar                | Per-site (gamma)                        | Per-partition relative rate |
| Branch length model       | Equal / proportional / unlinked                | Linked / scaled (default) / unlinked     | Single tree (Bayesian)                  | Linked or unlinked          |
| Partition count selection | ModelFinder greedy merge (BIC)                 | User-specified                           | Not applicable (infinite mixture)       | User-specified              |
| Convergence criterion     | BIC for model selection; ML for optimization   | Log-likelihood epsilon                   | bpcomp `maxdiff` + tracecomp `effsize`  | ESS in Tracer               |
| Memory layout             | Per-site x classes (C60: 112 GB); PMSF: 2.2 GB | Compressed patterns + per-partition CLVs | MCMC state: site allocations + profiles | XML model graph + BEAGLE    |
| Cost vs single-model      | C60: ~60x RAM; PMSF: ~1.2x                     | Scaled: ~1x per partition                | 10K-30K MCMC cycles                     | Scales with chain length    |

### 5.2 IQ-TREE PMSF

<a id="cite-8b"></a>[Wang et al. 2018](https://doi.org/10.1093/sysbio/syx068) [[8](#ref-8)] decouple profile estimation from tree search in a two-stage workflow:

1. Fit mixture model (e.g. LG+C60+F+G) on a guide tree with `-ft`, compute posterior mean site frequency profiles, write to `.sitefreq` file
2. Use fixed per-site profiles for tree search with `-fs`

Phase 1 can run on a high-memory machine with `-n 0` (stops after writing `.sitefreq`). Phase 2 requires only ~1.2x single-model RAM. PMSF enables nonparametric bootstrap under site-heterogeneous models.

IQ-TREE supports three partition modes: edge-equal (`-q`), edge-proportional (`-p`, recommended), edge-unlinked (`-Q`). ModelFinder's `MFP+MERGE` performs greedy partition merging using BIC.

### 5.3 RAxML-NG

Defaults to scaled branch lengths (shared topology + per-partition rate scalers). Three modes via `--brlen`: `linked`, `scaled` (default), `unlinked`. Per-partition model modifiers (`+G`, `+I`, `+R`) are independent. Per-rate-category likelihood scaling addresses numerical underflow for large trees ($> 5000$ taxa).

### 5.4 PhyloBayes CAT

Dirichlet process infinite mixture: the number of occupied frequency-profile components adapts to data. All components share the same exchangeability (CAT-Poisson uses uniform, CAT-GTR estimates from data). CAT-GTR has "virtually always the highest fit" except for very small datasets ($< 400$ positions). Best for 1K-20K positions, up to ~100 taxa.

Convergence requires 2+ independent chains: `maxdiff < 0.1` (bpcomp, tree topology) and `effsize > 300` (tracecomp, model parameters). `readpb_mpi -ss` exports posterior mean site-specific frequency profiles.

### 5.5 BEAST2

Independent link/unlink controls for substitution model, clock model, and tree model per partition. Default: unlinked substitution models, linked tree and clock. SiteModel wraps the substitution model with gamma rate categories. Computational cost scales linearly with category count. Package system allows novel models as plugins.

---

## 6. Model selection for partitioned analyses

### 6.1 Information criteria

<a id="cite-11"></a>[Posada and Buckley 2004](https://doi.org/10.1080/10635150490522304) [[11](#ref-11)] established AIC and Bayesian methods as preferred over hierarchical likelihood ratio tests. <a id="cite-12"></a>[Liu 2022](https://doi.org/10.1093/sysbio/syac081) [[12](#ref-12)] compared AIC and BIC on partition vs mixture models: AIC prefers complex mixture models and estimates branch lengths better; BIC prefers simpler models and estimates substitution parameters better. Cross-validation recommended as alternative.

### 6.2 Partition scheme selection

<a id="cite-13"></a>[Lanfear et al. 2012](https://doi.org/10.1093/molbev/mss020) [[13](#ref-13)] introduced PartitionFinder (greedy search with BIC/AIC). <a id="cite-14"></a>[Lanfear et al. 2014](https://doi.org/10.1186/1471-2148-14-82) [[14](#ref-14)] added relaxed hierarchical clustering for genome-scale datasets. <a id="cite-15"></a>[Frandsen et al. 2015](https://doi.org/10.1186/s12862-015-0283-7) [[15](#ref-15)] introduced k-means clustering of site rates as an automatic partitioning strategy, directly relevant to rate-based partitioning.

### 6.3 Implications for TreeTime

BIC is the safer default for partition count selection (penalizes overfitting, important for data-limited pathogen datasets). Cross-validation provides a model-free alternative but is computationally expensive.

---

## 7. Proposed approach for TreeTime

### 7.1 Partitioning strategies

Two complementary methods for assigning positions to partitions:

- Rate-based: run <a id="gloss-use-10"></a>Fitch <sup>[10](#gloss-10)</sup> parsimony reconstruction on the full alignment, count mutations per site, sort by mutation count, split into $q$ quantile bins. Related to the k-means site-rate clustering of <a id="cite-15b"></a>[Frandsen et al. 2015](https://doi.org/10.1186/s12862-015-0283-7) [[15](#ref-15)]
- Codon-position: parse GFF annotations to identify CDS regions, split protein-coding positions by 1st, 2nd, 3rd codon position, non-coding as 4th partition

<a id="gtr-family-inference"></a>

### 7.2 GTR family inference

A grouping specification controls which partitions share which parameters. The three tiers ($W$, $\pi$, $\mu$) form nested groups:

- $W$ groups: pool $n_{ij}$ and $\tau_j$ across partitions, update shared exchangeability
- $\pi$ groups: pool row-summed counts and root state frequencies, update shared equilibrium frequencies
- $\mu$ groups: pool total counts, update shared overall rate

### 7.3 Site-specific vs partitioned: distinct features

Site-specific models and discrete partitioning are related but not interchangeable. Site-specific GTR belongs as one dense partition with a `SiteSpecific` model variant, not as thousands of separate partitions. Representing per-site models as many partitions would lose $W$ sharing and destroy traversal efficiency. Discrete partitioning groups sites into coarser bins sharing all GTR parameters within a bin.

### 7.4 Molecular clock integration

The recommended workflow:

1. Fit initial scalar GTR and branch lengths
2. Run marginal reconstruction on dense profiles
3. Accumulate per-site $n_{ij}^a$, $\tau_j^a$, and root state
4. Update site-specific model via iterative inference
5. Normalize average rate to 1 (separates site model rate scale from molecular clock rate)
6. Re-run reconstruction and branch optimization with guarded convergence checks
7. Run clock inference with corrected branch lengths

Rate normalization at step 5 is essential: without it, the site model rate scale and molecular clock rate become confounded.

### 7.5 Two-stage estimation

Following IQ-TREE's PMSF pattern: estimate per-site preferences once from a good tree, then fix the preferences for subsequent analyses (rerooting, molecular clock, bootstrapping). This decouples expensive profile estimation from repeated tree operations. IQ-TREE serializes per-site frequencies to a `.sitefreq` text file (one line per site, whitespace-delimited frequency values) which can be transferred between machines and sessions. TreeTime would use a JSON format storing the full model ($W$, $\pi$, $\mu$, metadata) as its primary serialization, with `.sitefreq`-compatible export for cross-tool use. See the [serialization format proposal](../proposals/model-serialization-site-specific.md) for format details and examples.

---

## 8. Design recommendations

Informed by the software comparison and Puller et al. 2020 findings:

- R1: Add `enum GtrModel { Scalar(GTR), SiteSpecific(GTRSiteSpecific) }` in partition types. Value-based dispatch, no vtable overhead
- R2: Site-specific GTR is dense-only. Reject sparse mode when `--site-specific-gtr` is set. Sparse fixed-state aggregation is incompatible with per-site transition matrices
- R3: Encapsulate `Array2` vs `Array3` dimensionality behind model-level methods. Callers should not need to know which model variant is active
- R4: Add `get_mutation_counts_site_specific_dense()` alongside `get_mutation_counts_dense()`, preserving per-site $n_{ij}^a$ and $\tau_j^a$
- R5: Export inferred $\pi^a$, $\mu^a$, and shared $W$ in a serialized format for reuse as fixed profiles in later analyses. Primary format: JSON with ndarray serde (stores full model with metadata). Secondary: IQ-TREE-compatible `.sitefreq` text export for cross-tool interoperability (stores $\pi$ only). See [serialization format proposal](../proposals/model-serialization-site-specific.md)
- R6: Normalize average rate before clock use. Keep site model rate scale separate from molecular clock rate
- R7: Report diagnostics: mean $\delta\pi$, likelihood improvement, mean/max $\mu^a$, low-information site count, branch length change after model update. Add refusal thresholds for datasets with insufficient data
- R8: BIC as default criterion for partition count selection
- R9: Adopt RAxML-NG's scaled partition default (shared topology + per-partition rate scalers) rather than fully unlinked branch lengths

### 8.1 Approaches to avoid

- CAT/C60 as default pathogen model (targets deep phylogenomics, excessive memory for TreeTime's niche)
- CAT-BP lineage heterogeneity in core path (outside TreeTime's fast dated-pathogen niche)
- Unconstrained joint estimation of branch lengths, clock rate, $\mu^a$, and $\pi^a$ simultaneously (identifiability problems per Gascuel and Steel 2020)
- Per-site independent $W_{ij}^a$ (insufficient data support in pathogen genomics)

---

## 9. v1 codebase readiness

### 9.1 Existing infrastructure

| Component                                                       | Location                                                                                                                                   |
| --------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| `struct MutationCounts { nij, Ti, root_state }`                 | [packages/treetime/src/gtr/infer_gtr/common.rs#L10-L46](../../packages/treetime/src/gtr/infer_gtr/common.rs#L10-L46)                       |
| `fn infer_gtr_impl()` single-partition solver                   | [packages/treetime/src/gtr/infer_gtr/common.rs#L102-L163](../../packages/treetime/src/gtr/infer_gtr/common.rs#L102-L163)                   |
| `struct GTRSiteSpecific` with per-site eigendecomposition       | [packages/treetime/src/gtr/gtr_site_specific.rs#L47-L67](../../packages/treetime/src/gtr/gtr_site_specific.rs#L47-L67)                     |
| `fn infer_gtr_site_specific_impl()`                             | [packages/treetime/src/gtr/infer_gtr/site_specific.rs#L77-L187](../../packages/treetime/src/gtr/infer_gtr/site_specific.rs#L77-L187)       |
| `struct MutationCountsSiteSpecific { n_ija, T_ia, root_state }` | [packages/treetime/src/gtr/infer_gtr/site_specific.rs#L8-L26](../../packages/treetime/src/gtr/infer_gtr/site_specific.rs#L8-L26)           |
| Multi-partition `Vec<Arc<RwLock<...>>>`                         | [packages/treetime/src/partition/timetree.rs](../../packages/treetime/src/partition/timetree.rs)             |
| `trait HasGtr` (per-partition GTR access)                       | [packages/treetime/src/partition/traits.rs#L20-L32](../../packages/treetime/src/partition/traits.rs#L20-L32) |
| `fn PartitionFitch::compress()`                                 | [packages/treetime/src/partition/fitch.rs#L34-L56](../../packages/treetime/src/partition/fitch.rs#L34-L56)   |
| `fn Sub.pos()` for per-site mutation position                   | [packages/treetime/src/seq/mutation.rs#L21](../../packages/treetime/src/seq/mutation.rs#L21)                                               |
| `ExpQtInterpolator` (61-point non-uniform grid)                 | [packages/treetime/src/gtr/gtr_site_specific.rs#L433-L468](../../packages/treetime/src/gtr/gtr_site_specific.rs#L433-L468)                 |

The mathematical core for site-specific GTR is complete and well-tested: 4 golden master tests (1e-10 tolerance vs v0), 15 property tests (column stochastic, semigroup, equilibrium convergence, stationary distribution, interpolation accuracy). No production callers exist.

### 9.2 New code needed

| Component                                   | Notes                                                                                   |
| ------------------------------------------- | --------------------------------------------------------------------------------------- |
| `enum GtrModel { Scalar, SiteSpecific }`    | Replace `pub gtr: GTR` at 7+ partition call sites                                       |
| `get_mutation_counts_site_specific_dense()` | Per-site $n_{ij}^a$ and $\tau_j^a$ from tree traversal (port v0 `treeanc.py:1558-1589`) |
| Iterative refinement loop                   | Alternate model inference + ancestral reconstruction (port v0 `treeanc.py:1634-1677`)   |
| Per-site mutation counting from Fitch subs  | Iterate edges' `fitch_subs()`, count per `sub.pos()`                                    |
| Rate-based bin assignment                   | Sort positions by count, quantile split                                                 |
| Position-to-partition map type              | `Vec<usize>` or newtype                                                                 |
| Subsequence extraction from alignment       | Split FASTA records by position subset                                                  |
| GTR family inference (pooled groups)        | Wrap existing solver, decompose $W$/$\pi$/$\mu$ updates                                 |
| Partition grouping specification            | Data structure for $W$/$\pi$/$\mu$ group membership                                     |
| GFF/annotation parsing                      | No annotation infrastructure exists                                                     |
| Model serialization output                  | Export $\pi^a$, $\mu^a$, $W$ for reuse                                                  |
| Convergence diagnostics                     | Likelihood improvement, $\delta\pi$, rate statistics                                    |
| CLI integration                             | New flags/config for site-specific and partition strategy                               |

### 9.3 Related entries

- [kb/decisions/partition-system-architecture.md](../decisions/partition-system-architecture.md) - multi-partition design
- [kb/proposals/config-file-multi-partition.md](../proposals/config-file-multi-partition.md) - user-specified partition config
- [kb/issues/N-gtr-site-specific-partition-integration.md](../issues/N-gtr-site-specific-partition-integration.md) - site-specific GTR integration
- [kb/features/gtr.md](../features/gtr.md) - site-specific models listed as not implemented

---

## 10. Open questions

- No published analysis of convergence rate for the iterative EM-style update rules. The NMF connection guarantees monotone increase but says nothing about convergence speed or basin of attraction size
- The interaction between site-heterogeneous models and relaxed molecular clock estimation (lineage-rate variation) is not well-studied
- The optimal transition from per-site estimation to discrete partitioning (when to bin sites vs estimate individually) lacks formal treatment
- Epistatic interactions between sites (preferences at one site depending on state at another) are acknowledged but not modeled in any site-specific GTR framework

## Glossary

1. <a id="gloss-1"></a> **GTR (General Time Reversible).** The most general time-reversible nucleotide substitution model, parameterized by 6 exchangeability parameters ($W_{ij}$) and 4 stationary frequencies ($\pi_i$), scaled by an overall rate $\mu$ (<a id="cite-16"></a>[Tavare 1986](https://archive.org/details/someprobabilisticandstatisticalproblemsintheanalysisofdnasequences) [[16](#ref-16)]). [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Gamma rate variation.** Among-site rate heterogeneity modeled as a discretized Gamma distribution with shape parameter $\alpha$. Smaller $\alpha$ implies greater rate variation ([Yang 1994](https://doi.org/10.1007/BF00178256) [[3](#ref-3)]). [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **Mutation-selection balance.** Equilibrium between diversifying mutation and purifying selection at a genomic site, determining the site's stationary distribution of allowed states. The theoretical basis for site-specific equilibrium frequencies ([Halpern and Bruno 1998](https://doi.org/10.1093/oxfordjournals.molbev.a025995) [[1](#ref-1)]). [↩](#gloss-use-3)
4. <a id="gloss-4"></a> **Exchangeability matrix.** The symmetric matrix $W_{ij}$ in the GTR model encoding relative rates of exchange between nucleotide pairs, independent of base composition. [↩](#gloss-use-4)
5. <a id="gloss-5"></a> **Equilibrium frequencies.** The stationary distribution $\pi_i$ of nucleotide (or amino acid) frequencies under a substitution model, satisfying detailed balance. [↩](#gloss-use-5)
6. <a id="gloss-6"></a> **FreeRate (+R).** Site rate heterogeneity model with freely estimated rate categories and weights, without assuming a parametric distribution. More flexible than Gamma but requires more parameters. [↩](#gloss-use-6)
7. <a id="gloss-7"></a> **Dirichlet process.** A Bayesian nonparametric prior over probability distributions, used in the CAT model to define an infinite mixture of site frequency profiles where the number of occupied classes adapts to the data ([Lartillot and Philippe 2004](https://doi.org/10.1093/molbev/msh112) [[6](#ref-6)]). [↩](#gloss-use-7)
8. <a id="gloss-8"></a> **PMSF (Posterior Mean Site Frequency).** Per-site amino acid frequency profile computed as the conditional mean over mixture model classes, given data and a guide tree. Used as a fixed approximation to full mixture models for ML tree search ([Wang et al. 2018](https://doi.org/10.1093/sysbio/syx068) [[8](#ref-8)]). [↩](#gloss-use-8)
9. <a id="gloss-9"></a> **Non-negative matrix factorization (NMF).** A family of algorithms that decompose a non-negative matrix into non-negative factors using multiplicative update rules that guarantee positivity and monotone objective improvement ([Lee and Seung 2001](https://proceedings.neurips.cc/paper/2000/hash/f9d1152547c0bde01830b7e8bd60024c-Abstract.html) [[9](#ref-9)]). [↩](#gloss-use-9)
10. <a id="gloss-10"></a> **Fitch parsimony.** Ancestral state reconstruction algorithm minimizing the total number of state changes on a phylogenetic tree. Uses a two-pass approach: bottom-up to assign state sets, top-down to resolve ambiguities. [↩](#gloss-use-10)

## References

1. <a id="ref-1"></a> Halpern, Aaron L., and William J. Bruno. 1998. "Evolutionary Distances for Protein-Coding Sequences: Modeling Site-Specific Residue Frequencies." _Molecular Biology and Evolution_ 15(7):910-917. https://doi.org/10.1093/oxfordjournals.molbev.a025995 [↩¹](#cite-1a) [↩²](#cite-1b) [↩³](#cite-1c)
2. <a id="ref-2"></a> Puller, Vadim, Pavel Sagulenko, and Richard A. Neher. 2020. "Efficient Inference, Potential, and Limitations of Site-Specific Substitution Models." _Virus Evolution_ 6(2):veaa066. https://doi.org/10.1093/ve/veaa066 [↩¹](#cite-2a) [↩²](#cite-2b) [↩³](#cite-2c) [↩⁴](#cite-2d) [↩⁵](#cite-2e)
3. <a id="ref-3"></a> Yang, Ziheng. 1994. "Estimating the Pattern of Nucleotide Substitution." _Journal of Molecular Evolution_ 39:105-111. https://doi.org/10.1007/BF00178256 [↩¹](#cite-3) [↩²](#cite-3b) [↩³](#cite-3c)
4. <a id="ref-4"></a> Hilton, Sarah K., and Jesse D. Bloom. 2018. "Modeling Site-Specific Amino-Acid Preferences Deepens Phylogenetic Estimates of Viral Sequence Divergence." _Virus Evolution_ 4(2):vey033. https://doi.org/10.1093/ve/vey033 [↩](#cite-4)
5. <a id="ref-5"></a> Bruno, William J. 1996. "Modeling Residue Usage in Aligned Protein Sequences via Maximum Likelihood." _Molecular Biology and Evolution_ 13(10):1368-1374. https://doi.org/10.1093/oxfordjournals.molbev.a025583 [↩](#cite-5)
6. <a id="ref-6"></a> Lartillot, Nicolas, and Herve Philippe. 2004. "A Bayesian Mixture Model for Across-Site Heterogeneities in the Amino-Acid Replacement Process." _Molecular Biology and Evolution_ 21(6):1095-1109. https://doi.org/10.1093/molbev/msh112 [↩](#cite-6)
7. <a id="ref-7"></a> Lartillot, Nicolas, Henner Brinkmann, and Herve Philippe. 2007. "Suppression of Long-Branch Attraction Artefacts in the Animal Phylogeny Using a Site-Heterogeneous Model." _BMC Evolutionary Biology_ 7(Suppl 1):S4. https://doi.org/10.1186/1471-2148-7-S1-S4 [↩](#cite-7)
8. <a id="ref-8"></a> Wang, Huai-Chun, Bui Quang Minh, Edward Susko, and Andrew J. Roger. 2018. "Modeling Site Heterogeneity with Posterior Mean Site Frequency Profiles Accelerates Accurate Phylogenomic Estimation." _Systematic Biology_ 67(2):216-235. https://doi.org/10.1093/sysbio/syx068 [↩¹](#cite-8a) [↩²](#cite-8b)
9. <a id="ref-9"></a> Lee, Daniel D., and H. Sebastian Seung. 2001. "Algorithms for Non-Negative Matrix Factorization." In _Advances in Neural Information Processing Systems 13_, 556-562. MIT Press. https://proceedings.neurips.cc/paper/2000/hash/f9d1152547c0bde01830b7e8bd60024c-Abstract.html [↩](#cite-9)
10. <a id="ref-10"></a> Gascuel, Olivier, and Mike Steel. 2020. "A Darwinian Uncertainty Principle." _Systematic Biology_ 69(3):521-529. https://doi.org/10.1093/sysbio/syz054 [↩](#cite-10)
11. <a id="ref-11"></a> Posada, David, and Thomas R. Buckley. 2004. "Model Selection and Model Averaging in Phylogenetics: Advantages of Akaike Information Criterion and Bayesian Approaches Over Likelihood Ratio Tests." _Systematic Biology_ 53(5):793-808. https://doi.org/10.1080/10635150490522304 [↩](#cite-11)
12. <a id="ref-12"></a> Liu, Qin. 2022. "Performance of Akaike Information Criterion and Bayesian Information Criterion in Selecting Partition Models and Mixture Models." _Systematic Biology_ 72(2):469-479. https://doi.org/10.1093/sysbio/syac081 [↩](#cite-12)
13. <a id="ref-13"></a> Lanfear, Robert, Brett Calcott, Simon Y. W. Ho, and Stephane Guindon. 2012. "PartitionFinder: Combined Selection of Partitioning Schemes and Substitution Models for Phylogenetic Analyses." _Molecular Biology and Evolution_ 29(6):1695-1701. https://doi.org/10.1093/molbev/mss020 [↩](#cite-13)
14. <a id="ref-14"></a> Lanfear, Robert, Brett Calcott, David Kainer, Christoph Mayer, and Alexandros Stamatakis. 2014. "Selecting Optimal Partitioning Schemes for Phylogenomic Datasets." _BMC Evolutionary Biology_ 14:82. https://doi.org/10.1186/1471-2148-14-82 [↩](#cite-14)
15. <a id="ref-15"></a> Frandsen, Paul B., Brett Calcott, Christoph Mayer, and Robert Lanfear. 2015. "Automatic Selection of Partitioning Schemes for Phylogenetic Analyses Using Iterative k-Means Clustering of Site Rates." _BMC Evolutionary Biology_ 15:13. https://doi.org/10.1186/s12862-015-0283-7 [↩¹](#cite-15) [↩²](#cite-15b)
16. <a id="ref-16"></a> Tavare, Simon. 1986. "Some Probabilistic and Statistical Problems in the Analysis of DNA Sequences." _Lectures on Mathematics in the Life Sciences_ 17:57-86. https://archive.org/details/someprobabilisticandstatisticalproblemsintheanalysisofdnasequences [↩](#cite-16)
