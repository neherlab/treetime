# Tree inference for an end-to-end TreeTime pipeline

TreeTime infers time-scaled phylogenies, ancestral sequences, and molecular clocks on a given tree <a id="gloss-use-1"></a>topology <sup>[1](#gloss-1)</sup>. In Nextstrain and augur workflows, an external tree builder (IQ-TREE by default) makes this topology from the alignment, and TreeTime runs after it. This report surveys the modern techniques and tools for the tree-building stage, audits how the Nextstrain pipelines use them, maps the techniques onto TreeTime v1 components, and analyzes which techniques TreeTime can adopt to go from an alignment to a time tree without external programs.

Sequence alignment, the stage before tree building, is outside the scope of this report and is covered in [feat-aln.md](feat-aln.md). The report gives research findings and options. It does not approve a divergence from v0 or a scientific design; [Open decisions](#open-decisions) lists the decisions that the team must make.

## Executive summary

- **The pipelines use an external maximum-likelihood (ML) builder and then repair its output.** `augur tree` runs IQ-TREE with a reduced search that equals IQ-TREE's `-fast` mode. Pipelines then collapse near-zero branches into <a id="gloss-use-2"></a>polytomies <sup>[2](#gloss-2)</sup> (`--polytomy`), reroot in workflow code, or run a parsimony-based reversion and homoplasy fix that uses TreeTime v0 (mpox). Full-tree influenza builds already force IQ-TREE 3 into its CMAPLE mode. TreeTime v0 calls IQ-TREE, FastTree, or RAxML when no tree is given; TreeTime v1 requires `--tree`. See [How Nextstrain pipelines build trees](#how-nextstrain-pipelines-build-trees)
- **Outbreak data has low divergence, and the data determines the topology only weakly.** In a 64,000-genome SARS-CoV-2 tree, more than 70% of branches carry less than one substitution; 100 independent ML searches give trees with a mean pairwise Robinson-Foulds distance of 0.78, and about three quarters of them are statistically indistinguishable; 82% of inferred ancestral nodes have an identical sampled child. Parsimony placement with parsimony rearrangement (UShER with matOptimize) gives trees with likelihoods equal to or higher than IQ-TREE 2, FastTree 2, and RAxML-NG at a small fraction of their cost, and the sparse likelihood method MAPLE is more accurate than all of them. Five of the nine dataset families in `data/` pass the CMAPLE low-divergence test; dengue, Lassa, the H3N2 HA set, and the SNP-only tuberculosis set do not. See [Properties of the data](#properties-of-the-data-that-control-method-choice)
- **TreeTime v1 already contains most of the machinery for low-divergence data.** It has Fitch compression, a sparse marginal likelihood with inside and outside messages on every edge, mutation lists on edges, UShER mutation-annotated tree input and output, parsimony-guided polytomy moves, root-to-tip rerooting, and time-aware polytomy resolution. It has no starting-topology construction, no global rearrangement search, no incremental message update, no minimum-cost parsimony recurrence at multifurcations, and no branch support. See [TreeTime v1 building blocks and gaps](#treetime-v1-building-blocks-and-gaps)
- **The two-step pipeline has a measurable cost for dating.** When an IQ-TREE topology rooted by TreeTime or LSD was held fixed, only 62% of clock, node-age, and demographic parameters stayed within ±25% of the joint Bayesian estimate, against 93% for topologies taken from the Bayesian posterior. In 10 of 15 datasets, the ML topology contained none of the root splits of the posterior. Delphy shows that a timed tree with explicit mutations on its branches makes Bayesian joint inference feasible for 100,000 sequences. This survey found no published ML tree search that uses tip dates during the topology search; TreeTime's time-aware polytomy resolution is the closest existing mechanism. See [Time-aware topology inference](#time-aware-topology-inference)
- **All production ML search engines have GPL or AGPL licenses.** TreeTime (MIT license) can run them only as separate processes or reimplement their published algorithms. UShER and matOptimize, Nextclade, and Delphy are MIT-licensed references for low-divergence tree building. Deep-learning tree builders do not support nucleotide outbreak data at TreeTime's scale. See [Software ecosystem and licenses](#software-ecosystem-and-licenses)

## Background

### Pipeline stages

```mermaid
flowchart LR
  A["<b>Sequences</b><br/><small>FASTA</small>"] --> B["<b>Alignment</b><br/><small>augur align, Nextclade</small>"]
  B --> C["<b>Masking and filtering</b><br/><small>augur mask, augur filter</small>"]
  C --> D["<b>Tree building</b><br/><small>augur tree: IQ-TREE, FastTree, RAxML</small>"]
  D --> E["<b>Refinement and dating</b><br/><small>augur refine: TreeTime</small>"]
  E --> F["<b>Ancestral states and traits</b><br/><small>augur ancestral, augur traits: TreeTime</small>"]
  F --> G["<b>Export</b><br/><small>augur export: Auspice JSON</small>"]
  classDef external fill:#7a6a8a,stroke:#7a6a8a,color:#ffffff
  classDef treetime fill:#4a6a8a,stroke:#4a6a8a,color:#ffffff
  classDef other fill:#6b7b5e,stroke:#6b7b5e,color:#ffffff
  class D external
  class E,F treetime
  class A,B,C,G other
```

The TreeTime paper <a id="cite-1a"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[1](#ref-1)] states that "the time tree inference and dating are typically faster than the estimation of the tree topology", and that "tree building software often randomly resolves these polytomies into a series of bifurcations" in an order that is often inconsistent with the sampling dates. TreeTime therefore collapses zero-length branches and resolves the resulting polytomies with temporal information. Nextstrain <a id="cite-2"></a>[Hadfield et al. 2018](https://doi.org/10.1093/bioinformatics/bty407) [[2](#ref-2)] chains the stages above through augur <a id="cite-3"></a>[Huddleston et al. 2021](https://doi.org/10.21105/joss.02906) [[3](#ref-3)]. Tree building (purple) is the only stage of the core analysis that no Nextstrain-maintained program implements.

### What a tree builder optimizes

The tree-building methods in this report differ by the objective they optimize over topologies $T$:

$$
\hat T_{\text{pars}} = \arg\min_T \sum_{e \in T} \lvert M_e \rvert,
\qquad
\hat T_{\text{ML}} = \arg\max_{T} \max_{b,\theta} \sum_{s=1}^{L} \log P(D_s \mid T, b, \theta)
$$

where:

- $e$ -- an edge of $T$
- $M_e$ -- the set of substitutions that a most-parsimonious reconstruction puts on edge $e$
- $b$ -- the vector of branch lengths in expected substitutions per site
- $\theta$ -- the substitution model parameters, for example GTR rates and equilibrium frequencies
- $L$ -- the alignment length
- $D_s$ -- the alignment column at site $s$

Distance methods optimize a tree-length criterion over a matrix of pairwise distances. Bayesian methods sample topologies, branch lengths, and node dates together from a posterior distribution. TreeTime itself optimizes node dates and branch lengths on a fixed topology, and changes the topology only locally (polytomy resolution and rerooting).

### Why the topology stage matters for dating

<a id="cite-4a"></a>[Fourment et al. 2026](https://doi.org/10.1093/sysbio/syag069) [[4](#ref-4)] compared joint Bayesian inference in BEAST X with analyses that hold the topology fixed, on 15 viral datasets. One of the fixed-topology strategies is the standard two-step approach: an IQ-TREE ML topology, rooted by LSD or TreeTime. The authors report:

- Substitution and site-model parameters are robust to the fixed topology
- Clock rates, node ages, and population sizes change systematically. 93% of these parameters stayed within ±25% of the joint estimate for topologies from the Bayesian posterior, and 62% for IQ-TREE topologies rooted by TreeTime or LSD. The mean relative error of population sizes was nearly six times larger for the TreeTime and LSD rootings
- In 10 of 15 datasets, the IQ-TREE topology contained none of the root splits that the posterior sampled, so no rerooting of that topology can match the posterior root
- The authors call for "faster, time-aware methods that simultaneously integrate topology and parameter estimation"

Two limits apply. The study uses the unconstrained Bayesian analysis as the reference, which is itself a model-based approximation. It uses TreeTime and LSD only to select the root, and estimates node dates with BEAST on the fixed tree, so it measures the effect of the fixed topology and root, not TreeTime's dating.

## How Nextstrain pipelines build trees

The audit below reads the pipeline repositories at their revisions of 2026-09-21 to 2026-10-02.

### `augur tree`

- **Defaults.** `augur tree` adds fixed default arguments for each builder [[src](https://github.com/nextstrain/augur/blob/0d287496eed3816f94d674f5a08201273f479487/augur/tree.py#L29-L46)]: IQ-TREE `--ninit 2 -n 2 --epsilon 0.05 -T AUTO --redo` (the parts of IQ-TREE's `-fast` option), FastTree `-nt -nosupport`, RAxML `-f d -m GTRCAT -c 25 -p 235813`. The default method is IQ-TREE and the default model is GTR; `--substitution-model auto` runs IQ-TREE's model selection [[src](https://github.com/nextstrain/augur/blob/0d287496eed3816f94d674f5a08201273f479487/augur/tree.py#L466-L472)]. Augur version 34.1.4 (2026-09-09)
- **Name handling.** IQ-TREE rewrites some characters in sequence names, so augur replaces them before the call and restores them after [[src](https://github.com/nextstrain/augur/blob/0d287496eed3816f94d674f5a08201273f479487/augur/tree.py#L248-L264)]
- **VCF input.** For VCF input, augur writes only the parsimony-informative sites to the FASTA file that the builder reads [[src](https://github.com/nextstrain/augur/blob/0d287496eed3816f94d674f5a08201273f479487/augur/tree.py#L355-L415)] and adds no <a id="gloss-use-3"></a>ascertainment bias <sup>[3](#gloss-3)</sup> correction. The IQ-TREE documentation states that a `+ASC` model "should be applied if the alignment does not contain constant sites" [[doc](https://iqtree.github.io/doc/Substitution-Models#ascertainment-bias-correction)]. Branch lengths from this path are therefore in units of informative sites

### Settings of individual pipelines

- **ncov** keeps IQ-TREE as the only method and uses `-ninit 10 -n 4` [[src](https://github.com/nextstrain/ncov/blob/3432c85760b6c8cd0f815a85167ebb46c8131abb/defaults/parameters.yaml#L118-L119)], with IQ-TREE 2.2.0.3 pinned in its environment [[src](https://github.com/nextstrain/ncov/blob/3432c85760b6c8cd0f815a85167ebb46c8131abb/workflow/envs/nextstrain.yaml#L9)]. `augur refine` then runs with `--stochastic-resolve` unless polytomies are kept [[src](https://github.com/nextstrain/ncov/blob/3432c85760b6c8cd0f815a85167ebb46c8131abb/workflow/snakemake_rules/main_workflow.smk#L801)]
- **seasonal-flu** adds `-czb`, which collapses near-zero branches [[src](https://github.com/nextstrain/seasonal-flu/blob/f34feaa4650fe3cfd1fa41b11574f9647a92065a/profiles/nextstrain-public.yaml#L16)]. The full-tree profile uses `--pathogen-force -ninit 2 -n 2 --epsilon 0.05 -czb`, which forces IQ-TREE 3 to use the CMAPLE algorithm [[src](https://github.com/nextstrain/seasonal-flu/blob/f34feaa4650fe3cfd1fa41b11574f9647a92065a/profiles/full-trees.yaml#L20-L22)]
- **mpox** runs IQ-TREE with masked sites [[src](https://github.com/nextstrain/mpox/blob/3e484eb891e1b223f26029121fc94d1bbd891a6b/phylogenetic/rules/construct_phylogeny.smk#L34-L40)] and then an optional `fix_tree` step [[src](https://github.com/nextstrain/mpox/blob/3e484eb891e1b223f26029121fc94d1bbd891a6b/phylogenetic/rules/construct_phylogeny.smk#L43-L70)]. The script runs TreeTime v0 branch-length optimization with short-branch pruning under JC69 [[src](https://github.com/nextstrain/mpox/blob/3e484eb891e1b223f26029121fc94d1bbd891a6b/phylogenetic/scripts/fix_tree.py#L25-L26)], then repeats up to five times: it moves a grandchild that reverts a mutation of its parent node up to the grandparent, and it groups siblings that share mutations under a new node [[src](https://github.com/nextstrain/mpox/blob/3e484eb891e1b223f26029121fc94d1bbd891a6b/phylogenetic/scripts/fix_tree.py#L47-L130)]. The mpox Nextclade-dataset workflow uses `--polytomy` with an optional constraint tree [[src](https://github.com/nextstrain/mpox/blob/3e484eb891e1b223f26029121fc94d1bbd891a6b/nextclade/Snakefile#L518-L522)], and has an alternative CMAPLE rule (`-m JC --search EXHAUSTIVE --out-mul-tree --make-consistent`) [[src](https://github.com/nextstrain/mpox/blob/3e484eb891e1b223f26029121fc94d1bbd891a6b/nextclade/Snakefile#L482-L500)]
- **zika** uses the augur defaults [[src](https://github.com/nextstrain/zika/blob/e63cb9cd103fe51bf7ab804ff02c913842df4e81/phylogenetic/rules/construct_phylogeny.smk#L37-L41)]
- **ebola** replaces the defaults with `--polytomy --ninit 100 --epsilon 0.01` [[src](https://github.com/nextstrain/ebola/blob/eadbca0b2bbcd471835170486af3828c99093f13/phylogenetic/defaults/config.yaml#L97-L102)] and reroots on an outgroup with BioPython inside the workflow [[src](https://github.com/nextstrain/ebola/blob/eadbca0b2bbcd471835170486af3828c99093f13/phylogenetic/rules/construct_phylogeny.smk#L87-L107)]

### IQ-TREE 3 options that the pipelines use

IQ-TREE 3 <a id="cite-5a"></a>[Wong et al. 2026](https://doi.org/10.1093/molbev/msag117) [[5](#ref-5)], version 3.1.4 (2026-09-10):

- `-czb` and `--polytomy` are the same option [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/tools.cpp#L1923-L1926)]. It collapses internal branches shorter than four times the minimum branch length [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/main/phyloanalysis.cpp#L4003-L4006)]. The minimum is $10^{-6}$, or $0.1/L$ for alignments with $L \ge 100{,}000$ sites [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/main/phyloanalysis.cpp#L5259-L5266)]
- `--pathogen` lets IQ-TREE choose CMAPLE <a id="cite-6a"></a>[Ly-Trong et al. 2024](https://doi.org/10.1093/molbev/msae134) [[6](#ref-6)], and `--pathogen-force` always uses CMAPLE [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/tools.cpp#L2735-L2742)]. `--sprta` also forces CMAPLE [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/main/phyloanalysis.cpp#L5412-L5422)]
- The choice uses `fn cmaple::isEffective()` [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/main/phyloanalysis.cpp#L5490-L5512)]. CMAPLE encodes each sequence as a list of entries against a majority-consensus reference [[src](https://github.com/iqtree/cmaple/blob/3d45b1ab68e2d68a2825bf17a531e22200578cd6/alignment/alignment.cpp#L479-L530)] and accepts the alignment when no sequence has more than $0.067 L$ entries and the mean is at most $0.02 L$ [[src](https://github.com/iqtree/cmaple/blob/3d45b1ab68e2d68a2825bf17a531e22200578cd6/maple/cmaple.cpp#L17-L64)] [[src](https://github.com/iqtree/cmaple/blob/3d45b1ab68e2d68a2825bf17a531e22200578cd6/utils/tools.h#L357-L358)]

### TreeTime v0 and v1

- **v0** calls external programs when no tree is given. `fn tree_inference()` tries IQ-TREE, FastTree, and RAxML in that order as subprocesses ([`packages/legacy/treetime/treetime/utils.py#L410-L452`](../../packages/legacy/treetime/treetime/utils.py#L410-L452)); the IQ-TREE call uses HKY with `-ninit 2 -n 2 -me 0.05` ([`packages/legacy/treetime/treetime/utils.py#L481-L515`](../../packages/legacy/treetime/treetime/utils.py#L481-L515)). `fn assure_tree()` triggers it ([`packages/legacy/treetime/treetime/wrappers.py#L13-L22`](../../packages/legacy/treetime/treetime/wrappers.py#L13-L22))
- **v1** requires a tree. The `timetree` arguments declare `--tree` as optional, and the conversion step returns the required-argument error when it is missing ([`packages/app-commands/src/commands/timetree/args.rs#L157-L159`](../../packages/app-commands/src/commands/timetree/args.rs#L157-L159))
- **Stale KB premise.** [kb/issues/H-timetree-tree-inference-unimplemented.md](../issues/H-timetree-tree-inference-unimplemented.md) states that v0 infers a tree "using neighbor-joining or other methods". The inspected v0 code only calls external programs

### What the audit shows

- The pipelines tune the builder for speed (2 to 10 starting trees) and accept the result of a short search
- The pipelines carry tree post-processing that TreeTime functions already implement or can own: collapse of near-zero branches, rerooting, reversion and homoplasy repair, and time-aware polytomy resolution
- The only use of a parsimony-type or sparse-likelihood builder is CMAPLE, through IQ-TREE 3 (full influenza trees) or directly (mpox Nextclade dataset, optional)

## Properties of the data that control method choice

### Low divergence and weak topological signal

For an edge with length $b_e$, the expected number of substitutions on the edge is $\lambda_e = L b_e$. Under a Poisson model, the probability that the edge carries no substitution at all is $e^{-\lambda_e}$. When $\lambda_e < 1$, which is common in outbreak data, most such edges carry zero or one substitution, and the alignment cannot resolve the branching order around them. Three measurements show the extent:

- In the 64,000-genome SARS-CoV-2 benchmark of <a id="cite-7a"></a>[Wang et al. 2023](https://doi.org/10.1093/bioinformatics/btad536) [[7](#ref-7)], more than 70% of branches represent less than one substitution, and "most are resolved effectively at random"
- <a id="cite-8a"></a>[Morel et al. 2021](https://doi.org/10.1093/molbev/msaa314) [[8](#ref-8)] ran 100 independent RAxML-NG searches on 4,869 SARS-CoV-2 genomes. The resulting trees had a mean pairwise relative Robinson-Foulds distance of 0.78, and 74 to 76 of the 100 trees (depending on the alignment version) passed the statistical tests for plausible trees. Searches from parsimony starting trees scored more than 400 log-likelihood units better than searches from random trees
- In a multifurcating SARS-CoV-2 tree of 364,427 samples, 68,261 of 83,216 inferred ancestral nodes (82%) have a sampled child with an identical genotype, that is, a zero-length branch <a id="cite-9a"></a>[Kramer et al. 2023](https://doi.org/10.1093/sysbio/syad031) [[9](#ref-9)]

The TreeTime paper reports that TreeTime "was tested predominantly on sequences from viruses with a pairwise identity above 90%" <a id="cite-1b"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[1](#ref-1)]. TreeTime's main use is therefore in the same regime.

### Bootstrap support of single-mutation branches

Consider a branch whose only evidence is $k$ sites that carry its mutation, in an alignment without conflicting characters. A bootstrap replicate draws $L$ sites with replacement. The replicate recovers the branch when it draws at least one of the $k$ sites:

$$
P(\text{branch recovered}) = 1 - \left(1 - \frac{k}{L}\right)^{L} \approx 1 - e^{-k}
$$

For $k = 1$ this gives $1 - e^{-1} \approx 0.63$, and for $k = 3$ it gives $1 - e^{-3} \approx 0.95$. The two values agree with the published figures: about 63% expected Felsenstein bootstrap support for a single-mutation branch <a id="cite-10a"></a>[Lemoine and Gascuel 2024](https://doi.org/10.1093/molbev/msae238) [[10](#ref-10)], and about three mutations for 95% support <a id="cite-11a"></a>[De Maio et al. 2025](https://doi.org/10.1038/s41586-025-09567-x) [[11](#ref-11)]. A correct branch with one supporting mutation thus looks weak under the standard bootstrap. See [Branch support](#branch-support) for alternatives.

### Divergence of the TreeTime example datasets

To estimate which datasets fit low-divergence methods, the following test applies the CMAPLE suitability rule to the largest alignment of each family in `data/`. The method:

- Reference: the per-site majority base among A, C, G, T (CMAPLE also builds a majority-consensus reference when none is given)
- Entries per sequence: substitutions against the reference, plus ambiguous bases, plus one entry for each run of `N` or gap characters, as in CMAPLE's encoding
- Gate: maximum entries per site at most 0.067 and mean entries per site at most 0.02
- The p-distance column counts substitutions over sites where both the sequence and the reference have a called base

| Dataset                        | Sequences |   Sites | Mean entries/site | Max entries/site | Mean p-distance | Max p-distance | Passes gate |
| ------------------------------ | --------: | ------: | ----------------: | ---------------: | --------------: | -------------: | ----------- |
| `flu/h3n2/500` (HA)            |       476 |   1,407 |            0.0222 |           0.1151 |          0.0222 |         0.1151 | no          |
| `ebola/362`                    |       362 |  19,006 |            0.0034 |           0.0194 |          0.0006 |         0.0019 | yes         |
| `zika/86`                      |        86 |  10,807 |            0.0020 |           0.0037 |          0.0020 |         0.0037 | yes         |
| `tb/149` (variable sites only) |       149 |     216 |            0.1279 |           0.1806 |          0.1275 |         0.1806 | no          |
| `rsv/a/2000`                   |     1,999 |  15,225 |            0.0136 |           0.0631 |          0.0136 |         0.0631 | yes         |
| `dengue/2000`                  |     1,954 |  10,723 |            0.0487 |           0.1753 |          0.0494 |         0.1811 | no          |
| `lassa/L/500`                  |       414 |  10,402 |            0.0994 |           0.2212 |          0.1213 |         0.3013 | no          |
| `mpox/clade-ii/2000`           |     1,989 | 197,209 |            0.0003 |           0.0101 |          0.0001 |         0.0035 | yes         |
| `sc2/4500`                     |     5,100 |  29,903 |            0.0024 |           0.0294 |          0.0020 |         0.0065 | yes         |

The result divides the datasets into two groups:

- **Low divergence** (Ebola, Zika, RSV-A, mpox clade II, SARS-CoV-2): the regime in which parsimony placement and sparse likelihood methods were developed and benchmarked
- **Divergent or special** (H3N2 HA over many seasons, dengue, Lassa): methods built on short-branch approximations lose accuracy or speed; ML search or distance methods are the established approach. The tuberculosis alignment contains only 216 variable sites. Its per-site divergence has no meaning without the number of invariant sites, and a likelihood on it needs constant-site counts or an ascertainment correction

These figures are a proxy for CMAPLE's own test. They were computed with a separate script that mirrors CMAPLE's counting rule, and small differences in the handling of ambiguous characters can change the values slightly.

## Landscape of tree-building methods

### Distance methods

- **Neighbor joining (NJ)** <a id="cite-12"></a>[Saitou and Nei 1987](https://doi.org/10.1093/oxfordjournals.molbev.a040454) [[12](#ref-12)] joins, at each step, the pair of nodes that minimizes the total tree length, starting from a star tree. The standard algorithm takes $O(n^3)$ time and $O(n^2)$ memory for $n$ sequences. NJ is a greedy algorithm for the <a id="gloss-use-4"></a>balanced minimum evolution <sup>[4](#gloss-4)</sup> (BME) criterion <a id="cite-13"></a>[Gascuel and Steel 2006](https://doi.org/10.1093/molbev/msl072) [[13](#ref-13)]
- **BIONJ** <a id="cite-14"></a>[Gascuel 1997](https://doi.org/10.1093/oxfordjournals.molbev.a025808) [[14](#ref-14)] uses a variance model of the distances in the reduction step at the same cost as NJ. Its gain over NJ is small at low divergence and large at high and variable rates
- **BME and FastME.** <a id="cite-15"></a>[Desper and Gascuel 2002](https://doi.org/10.1089/106652702761034136) [[15](#ref-15)] build a BME tree in $O(n^2 \cdot \mathrm{diam}(T))$ time and improve it with nearest-neighbor interchanges, where $\mathrm{diam}(T)$ is the topological diameter of the tree. FastME 2.0 adds subtree pruning and regrafting moves and stays as fast as NJ <a id="cite-16"></a>[Lefort, Desper, and Gascuel 2015](https://doi.org/10.1093/molbev/msv150) [[16](#ref-16)]
- **Faster exact NJ.** RapidNJ <a id="cite-17"></a>[Simonsen, Mailund, and Pedersen 2008](https://doi.org/10.1007/978-3-540-87361-7_10) [[17](#ref-17)] prunes the search for the pair to join. DecentTree <a id="cite-7b"></a>[Wang et al. 2023](https://doi.org/10.1093/bioinformatics/btad536) [[7](#ref-7)] vectorizes NJ and BIONJ; on 64,000 SARS-CoV-2 sequences it was 1.8 times (1 thread) and 5.6 times (32 threads) faster than RapidNJ, and needed about 50 GB of memory against 17 GB for FastTree. The authors state that it "may not be applicable" to millions of sequences
- **NJ at larger scale.** Dynamic and heuristic NJ reach one million taxa <a id="cite-18"></a>[Clausen 2023](https://doi.org/10.1093/bioinformatics/btac774) [[18](#ref-18)]. Sparse NJ computes $O(n \log n)$ distance entries instead of all $n^2$ <a id="cite-19"></a>[Kurt, Bouchard-Côté, and Lagergren 2024](https://doi.org/10.1093/bioinformatics/btae701) [[19](#ref-19)]. FastTree's profile-based NJ avoids the distance matrix and uses $O(NLa)$ memory for $N$ sequences, $L$ sites, and alphabet size $a$ <a id="cite-20"></a>[Price, Dehal, and Arkin 2009](https://doi.org/10.1093/molbev/msp077) [[20](#ref-20)]
- **Multifurcating NJ** <a id="cite-21"></a>[Fernández, Segura-Alabart, and Serratosa 2023](https://doi.org/10.1007/s00239-023-10134-z) [[21](#ref-21)] outputs polytomies where distances do not separate the candidates, so the result does not depend on the input order
- **Distances.** Corrected distances (JC69, K2P, TN93 <a id="cite-22"></a>[Tamura and Nei 1993](https://doi.org/10.1093/oxfordjournals.molbev.a040023) [[22](#ref-22)], or ML distances) need explicit rules for gaps, ambiguous bases, unequal coverage, saturation, and zero distances. IQ-TREE skips a site in a pair when either sequence has a gap or an ambiguous character there, and assigns the maximum distance 9.0 to pairs without overlap or with saturation [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/tools.h#L531)]

**Fit for TreeTime.** NJ or BIONJ gives a complete, deterministic starting tree from any alignment, including divergent ones. Its limits are the $O(n^2)$ distance matrix and its behavior at low divergence, where it resolves short branches arbitrarily and returns a binary tree.

### Parsimony and placement on mutation-annotated trees

- **Fitch and Sankoff.** Fitch's algorithm <a id="cite-23"></a>[Fitch 1971](https://doi.org/10.1093/sysbio/20.4.406) [[23](#ref-23)] finds minimum-change ancestral states on a fixed binary tree; Sankoff's algorithm <a id="cite-24"></a>[Sankoff 1975](https://doi.org/10.1137/0128004) [[24](#ref-24)] generalizes it to any cost matrix. The intersection-or-union recurrence is exact only for binary nodes; the recurrence for nodes of any degree keeps the states with minimum subtree cost <a id="cite-25"></a>[Hartigan 1973](https://doi.org/10.2307/2529676) [[25](#ref-25)]
- **UShER** <a id="cite-26a"></a>[Turakhia et al. 2021](https://doi.org/10.1038/s41588-021-00862-7) [[26](#ref-26)] stores a <a id="gloss-use-5"></a>mutation-annotated tree <sup>[5](#gloss-5)</sup> (MAT) and places each new sample at the node where the extra parsimony cost is smallest. The cost at a node is the difference between the sample's mutations and the mutations on the path from the root. With a preprocessed MAT, one sample takes about 0.5 s, against about 28 CPU minutes and 791 GB for EPA-ng. On simulated data, UShER finds the correct sister node for 97.2% of samples (98.5% when the placement is unique). Ties and the number of equally parsimonious placements are reported
- **matOptimize** <a id="cite-27"></a>[Ye et al. 2022](https://doi.org/10.1093/bioinformatics/btac401) [[27](#ref-27)] improves a MAT with parallel parsimony <a id="gloss-use-6"></a>subtree pruning and regrafting <sup>[6](#gloss-6)</sup> (SPR) moves. It computes the score change of a move incrementally from Fitch sets near the move, and doubles the SPR radius in each round. It optimized a 3-million-sample tree in 7.7 h. The paper states a default stopping threshold of 0.5% improvement per round; the code default is 0.0005 (0.05%) [[src](https://github.com/yatisht/usher/blob/ac9c982d937c3bc43e2c16a0b73d30cf2b937118/src/matOptimize/main.cpp#L169)]
- **Online parsimony against ML.** <a id="cite-9b"></a>[Kramer et al. 2023](https://doi.org/10.1093/sysbio/syad031) [[9](#ref-9)] compared online UShER with matOptimize to IQ-TREE 2, FastTree 2, and RAxML-NG on SARS-CoV-2. At 26,486 samples, the parsimony scores were 16,130 (matOptimize), 16,195 (IQ-TREE 2), and 16,290 (FastTree 2); matOptimize needed 6 s and 0.15 GB, IQ-TREE 2 needed 3 h 29 min and 72 GB. With up to 14 days of runtime on 4,500 to 13,200 samples, UShER with matOptimize gave higher log-likelihoods than IQ-TREE 2 and RAxML-NG on real data
- **Uncertainty and order effects.** matUtils reports the number of equally parsimonious placements and the neighborhood size of each sample <a id="cite-28a"></a>[McBroome et al. 2021](https://doi.org/10.1093/molbev/msab264) [[28](#ref-28)]. For the 4.47-million-sample Viridian tree, the order in which samples enter UShER changed the deep structure of the tree; random order misplaced variants of concern, and the authors first added complete samples in temporal order <a id="cite-29"></a>[Hunt et al. 2026](https://doi.org/10.1038/s41592-025-02947-1) [[29](#ref-29)]
- **Nextclade** <a id="cite-30"></a>[Aksamentov et al. 2021](https://doi.org/10.21105/joss.03773) [[30](#ref-30)] places query sequences on a reference tree by scanning all nodes for the smallest mutation distance [[src](https://github.com/nextstrain/nextclade/blob/529db53e8f5b16d8a5d1161abd8b60ef04109817/packages/nextclade/src/tree/tree_find_nearest_node.rs#L21)]. Its optional greedy tree builder adds queries in order of increasing private mutation count [[src](https://github.com/nextstrain/nextclade/blob/529db53e8f5b16d8a5d1161abd8b60ef04109817/packages/nextclade/src/tree/tree_builder.rs#L29-L32)], moves each query to a neighboring node while that reduces private mutations, and splits branches at shared mutations. Its documentation states that this placement does not replace a full phylogenetic analysis [[doc](https://docs.nextstrain.org/projects/nextclade/en/stable/user/algorithm/03-phylogenetic-placement.html)]
- **Delphy's tree initializer** <a id="cite-31a"></a>[Varilly et al. 2026](https://doi.org/10.1038/s41586-026-11012-6) [[31](#ref-31)] builds a complete starting tree from a reference sequence and the tip differences, without external programs [[src](https://github.com/broadinstitute/delphy/blob/5d5989f3e4cac9ed5ae9bbcad62bcbc81afd8c2a/core/utree.h#L232-L326)]: a guide tree that inserts each tip at the edge with the fewest new differences (branch-and-bound search), a refined tree that adds tips in nearest-first order of the guide tree, parsimony SPR refinement, and rooting at the position that maximizes the root-to-tip regression fit against tip dates (OLS or GLS). Delphy is MIT-licensed and runs in the browser
- **Limits of parsimony.** Parsimony is statistically inconsistent when parallel changes outnumber informative changes (<a id="gloss-use-7"></a>long-branch attraction <sup>[7](#gloss-7)</sup>) <a id="cite-32"></a>[Felsenstein 1978](https://doi.org/10.1093/sysbio/27.4.401) [[32](#ref-32)]. In the MAPLE benchmarks <a id="cite-33a"></a>[De Maio et al. 2023](https://doi.org/10.1038/s41588-023-01368-0) [[33](#ref-33)], matOptimize was less accurate than ML methods on simulated data, more accurate on real data, second only to MAPLE, and its relative performance decreased at higher divergence

**Fit for TreeTime.** Parsimony placement and parsimony SPR operate on the same data that TreeTime's sparse partitions hold: the root sequence and the substitutions on each edge. The evidence supports them in the low-divergence group of datasets and does not support them for the divergent group.

### Maximum-likelihood search

- **IQ-TREE.** IQ-TREE combines hill climbing with stochastic perturbation <a id="cite-34"></a>[Nguyen et al. 2015](https://doi.org/10.1093/molbev/msu300) [[34](#ref-34)]. The v3.1.4 defaults build 100 parsimony starting trees by randomized stepwise addition, run <a id="gloss-use-8"></a>nearest-neighbor interchange <sup>[8](#gloss-8)</sup> (NNI) hill climbing on the best 20, and keep a candidate set of 5 [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/tools.cpp#L7421-L7434)] [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/tools.cpp#L7227)] [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/tree/iqtree.cpp#L836-L839)]. `-fast` reduces this to 2 starting trees and 2 iterations [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/tools.cpp#L4196-L4213)]. IQ-TREE removes identical sequences before the search by default [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/tools.cpp#L7468)]. IQ-TREE 2 extended the model set and the parallel computation for genomic data <a id="cite-35"></a>[Minh et al. 2020](https://doi.org/10.1093/molbev/msaa015) [[35](#ref-35)]. IQ-TREE 3 <a id="cite-5b"></a>[Wong et al. 2026](https://doi.org/10.1093/molbev/msag117) [[5](#ref-5)] adds mixture models, concordance factors, integration with divergence time estimation, and a sequence simulator; the source also integrates CMAPLE and SPRTA (see [IQ-TREE 3 options that the pipelines use](#iq-tree-3-options-that-the-pipelines-use))
- **RAxML-NG** <a id="cite-36"></a>[Kozlov et al. 2019](https://doi.org/10.1093/bioinformatics/btz305) [[36](#ref-36)] reimplements the greedy SPR search of RAxML. Adaptive RAxML-NG <a id="cite-37"></a>[Togkousidis et al. 2023](https://doi.org/10.1093/molbev/msad227) [[37](#ref-37)] sets the search effort from a predicted dataset difficulty (Pythia <a id="cite-38a"></a>[Haag et al. 2022](https://doi.org/10.1093/molbev/msac254) [[38](#ref-38)]). RAxML-NG 2.0 (2.0.3, 2026-09-03) makes the adaptive search the default and adds a fast mode with early stopping; with machine-learning branch support prediction, the fast mode reduces inference time 65-fold against RAxML-NG 1.2 <a id="cite-39"></a>[Kozlov et al. 2026](https://doi.org/10.64898/2026.09.09.750097) [[39](#ref-39)]. The minimum branch length is $10^{-6}$ [[src](https://github.com/amkozlov/raxml-ng/blob/d396351ee1263b704adb135edd6a4a84522dbd5b/src/constants.hpp#L29)]; when the best tree contains near-zero branches, RAxML-NG writes a warning and an additional tree with these branches collapsed [[src](https://github.com/amkozlov/raxml-ng/blob/d396351ee1263b704adb135edd6a4a84522dbd5b/src/main.cpp#L3026-L3043)]
- **FastTree 2** <a id="cite-40"></a>[Price, Dehal, and Arkin 2010](https://doi.org/10.1371/journal.pone.0009490) [[40](#ref-40)] starts from profile NJ, improves the tree with minimum-evolution NNI and SPR moves, then with ML NNI moves, and approximates rate variation with one rate category per site (CAT). It is 100 to 1,000 times faster than full ML SPR search on large alignments, with most disagreeing splits poorly supported. VeryFastTree 4 parallelizes the same algorithm and builds a one-million-sequence tree in 36 h on one server <a id="cite-41"></a>[Piñeiro and Pichel 2024](https://doi.org/10.1093/gigascience/giae055) [[41](#ref-41)]
- **PhyML 3** <a id="cite-42"></a>[Guindon et al. 2010](https://doi.org/10.1093/sysbio/syq010) [[42](#ref-42)] adds an SPR search of user-defined intensity, in which the parsimony criterion filters out the least promising moves
- **Benchmarks on general data.** On single-gene alignments of 19 phylogenomic datasets, IQ-TREE and RAxML found the highest likelihood in 80.17% and 75.99% of alignments and FastTree in 1.67%, with FastTree much faster <a id="cite-43"></a>[Zhou et al. 2018](https://doi.org/10.1093/molbev/msx302) [[43](#ref-43)]
- **Pathologies on low-divergence data.** The likelihood surface is flat with many optima (see [Low divergence and weak topological signal](#low-divergence-and-weak-topological-signal)). NNI search can stop early in multifurcating regions of tree space, which dominate real data, while SPR moves escape them <a id="cite-44"></a>[Whelan and Money 2010](https://doi.org/10.1093/molbev/msq163) [[44](#ref-44)]. All ML builders return binary trees with branch lengths clamped at a minimum value, so the true polytomies of the data appear as chains of near-zero branches in arbitrary order

**Model selection.** ModelFinder <a id="cite-45"></a>[Kalyaanamoorthy et al. 2017](https://doi.org/10.1038/nmeth.4285) [[45](#ref-45)] fits candidate substitution and rate models and selects one by an information criterion. Topologies inferred under different selected models, or simply under GTR with rate variation, are usually very similar <a id="cite-46"></a>[Abadi et al. 2019](https://doi.org/10.1038/s41467-019-08822-w) [[46](#ref-46)]. On SARS-CoV-2 data, the FreeRate model gave unstable likelihood scores, and <a id="cite-8b"></a>[Morel et al. 2021](https://doi.org/10.1093/molbev/msaa314) [[8](#ref-8)] recommend the discrete gamma model for such data.

**Fit for TreeTime.** A full ML search engine is the established method for divergent data. The production engines are large GPL or AGPL code bases (see [Software ecosystem and licenses](#software-ecosystem-and-licenses)), and their search heuristics assume dense likelihood evaluation. Their reported accuracy is the reference that an internal builder is measured against.

### Pandemic-scale likelihood: MAPLE and CMAPLE

MAPLE <a id="cite-33b"></a>[De Maio et al. 2023](https://doi.org/10.1038/s41588-023-01368-0) [[33](#ref-33)] replaces partial likelihood vectors with <a id="gloss-use-9"></a>genome lists <sup>[9](#gloss-9)</sup>: run-length entries relative to a reference genome, where one entry covers a run of reference positions and only differing or uncertain positions get explicit vectors. It also replaces the matrix exponential with a first-order approximation:

$$
P(t) = e^{tQ} \approx I + tQ
$$

where $Q$ is the instantaneous rate matrix and $t$ the branch length. MAPLE warns when an estimated branch length exceeds 0.01 or a genome differs from the reference by more than 10%. It builds the tree by stepwise placement and then improves it with SPR moves, updating genome lists only in a local part of the tree after each change. Reported results:

- Higher accuracy than RAxML-NG, the most accurate of the compared methods, on simulated and real SARS-CoV-2 data, while being more than 100-fold faster
- Trees about 25 times larger than IQ-TREE 2 or FastTree 2 can handle (500,000 against 20,000 samples); 500,000 samples in 69.4 h and 8.4 GB on one core
- At about 50 times the divergence of SARS-CoV-2, traditional ML methods become more efficient than MAPLE, while MAPLE's accuracy stays high

**Branch lengths under the first-order approximation.** With the approximation, the likelihood of site $i$ on a branch of length $t$ is linear in $t$: $\mathcal{L}_i(t) = u_i^{\top}(I + tQ)\,v_i = \alpha_i + \beta_i t$, where $u_i$ and $v_i$ are the partial likelihood vectors above and below the branch, $\alpha_i = u_i^{\top} v_i$, and $\beta_i = u_i^{\top} Q v_i$. Sites where both normalized vectors concentrate on the same state ($\alpha_i = 1$) contribute $\log(1 + \beta_i t) \approx \beta_i t$ for small $t$, and sites where they differ contribute $\log \beta_i + \log(c_i + t)$ with $c_i = \alpha_i / \beta_i$. The branch log-likelihood therefore has the form

$$
\ell(t) = a\,t + \sum_{i \in \mathcal{D}} \log(c_i + t) + \text{const},
\qquad
\ell'(t) = a + \sum_{i \in \mathcal{D}} \frac{1}{c_i + t}
$$

where $a$ collects the linear terms of the agreeing sites and $\mathcal{D}$ is the set of differing sites. CMAPLE finds the root of $\ell'(t)$ numerically in this form [[src](https://github.com/iqtree/cmaple/blob/3d45b1ab68e2d68a2825bf17a531e22200578cd6/tree/tree.cpp#L7603-L7720)].

**CMAPLE and later work.**

- CMAPLE <a id="cite-6b"></a>[Ly-Trong et al. 2024](https://doi.org/10.1093/molbev/msae134) [[6](#ref-6)] reimplements MAPLE in C++. On 10,000 SARS-CoV-2 sequences it took 8 min and 0.24 GB, against 3 h for FastTree 2 and 40.5 h for IQ-TREE 2, with higher likelihoods. The tree code alone has about 15,000 lines (`tree/tree.cpp`). Version 2.0.0 (2026-05-01), GPL-2.0
- Rate variation and sequencing errors: MAPLE adds one rate per site and a per-site error probability, estimated by expectation maximization, and uses tree nodes as local references so that divergence from a single global reference does not grow with tree depth <a id="cite-47"></a>[De Maio et al. 2026](https://doi.org/10.1038/s41592-025-02932-8) [[47](#ref-47)]
- CMAPLE 2 (preprint) adds multithreaded placement, SPR search, and SPRTA, multiple references, and MAT output <a id="cite-48"></a>[Ly-Trong et al. 2026](https://doi.org/10.64898/2026.06.15.732229) [[48](#ref-48)]

**Comparison with TreeTime's sparse likelihood.** TreeTime's sparse marginal design names MAPLE as the model and notes that a MAPLE genome list stops being sparse as the tree gets deeper ([kb/algo/ancestral.md](../algo/ancestral.md), "Sparse marginal design"). The v1 implementation differs from MAPLE in three ways:

- It keeps an explicit vector only for positions that are variable in the local subtree, and represents all other positions by one vector per state with a count of positions (`struct SparseSeqDistribution` [`packages/treetime/src/partition/storage/sparse.rs#L155`](../../packages/treetime/src/partition/storage/sparse.rs#L155)). The variable set follows the local subtree, so it does not grow with the distance to a global reference
- It propagates messages with the exact matrix exponential (`exp_qt.dot(...)` in [`packages/treetime/src/partition/marginal/sparse/message.rs#L93-L130`](../../packages/treetime/src/partition/marginal/sparse/message.rs#L93-L130)), with no short-branch approximation
- It keeps the contribution of positions without variation: `fn combine_messages()` adds the count of each fixed state times its log normalization to the log-likelihood ([`packages/treetime/src/partition/marginal/sparse/message.rs#L88`](../../packages/treetime/src/partition/marginal/sparse/message.rs#L88)). An ascertainment correction is therefore needed only for inputs that contain no invariant sites, such as the variable-site tuberculosis alignment or augur's VCF path

The v1 sparse passes drop a variable position from a message when all child messages agree and the peak probability is above a threshold ([`packages/treetime/src/partition/marginal/sparse/message.rs#L21`](../../packages/treetime/src/partition/marginal/sparse/message.rs#L21)). This is an approximation, and [kb/issues/M-ancestral-dense-sparse-divergence.md](../issues/M-ancestral-dense-sparse-divergence.md) reports dense-sparse log-likelihood differences of up to $7.6 \times 10^{-6}$ relative on about 2.5% of random configurations.

### Placement on a fixed backbone tree

<a id="gloss-use-10"></a>Phylogenetic placement <sup>[10](#gloss-10)</sup> adds query sequences to an existing reference tree without changing it:

- **pplacer** <a id="cite-49"></a>[Matsen, Kodner, and Armbrust 2010](https://doi.org/10.1186/1471-2105-11-538) [[49](#ref-49)] precomputes the partial likelihoods on both sides of every edge in two traversals, then evaluates each candidate edge as a three-taxon tree, in time linear in the number of taxa and the query length
- **EPA-ng** <a id="cite-50"></a>[Barbera et al. 2019](https://doi.org/10.1093/sysbio/syy054) [[50](#ref-50)] uses a fast preselection of candidate edges followed by a thorough evaluation, and scales to billions of short reads on a fixed tree
- **APPLES** <a id="cite-51"></a>[Balaban, Sarmashghi, and Mirarab 2020](https://doi.org/10.1093/sysbio/syz063) [[51](#ref-51)] places queries from distances and supports backbones of about 200,000 leaves
- **SCAMPP** <a id="cite-52"></a>[Wedell, Cai, and Warnow 2023](https://doi.org/10.1109/tcbb.2022.3170386) [[52](#ref-52)] restricts each placement to a query-specific subtree so that likelihood placement scales to large backbones

For SARS-CoV-2 samples on a 38,342-leaf tree, UShER placed one sample in about 0.5 s and EPA-ng needed about 28 CPU minutes and 791 GB of memory <a id="cite-26b"></a>[Turakhia et al. 2021](https://doi.org/10.1038/s41588-021-00862-7) [[26](#ref-26)]. Placement alone needs a backbone tree, so it cannot build the first tree; it fits incremental updates of a saved tree.

### Branch support

- **Felsenstein bootstrap** <a id="cite-53"></a>[Felsenstein 1985](https://doi.org/10.1111/j.1558-5646.1985.tb00420.x) [[53](#ref-53)] repeats the search on resampled alignments. It costs one full search per replicate and underestimates the support of short correct branches (see [Bootstrap support of single-mutation branches](#bootstrap-support-of-single-mutation-branches))
- **UFBoot2** <a id="cite-54"></a>[Hoang et al. 2018](https://doi.org/10.1093/molbev/msx281) [[54](#ref-54)] approximates the bootstrap during the search and corrects its overestimation of support on polytomies
- **SH-aLRT and aBayes** test each branch locally against its NNI alternatives <a id="cite-55"></a>[Anisimova and Gascuel 2006](https://doi.org/10.1080/10635150600755453) [[55](#ref-55)]; <a id="cite-56"></a>[Anisimova et al. 2011](https://doi.org/10.1093/sysbio/syr041) [[56](#ref-56)]. Their values are not bootstrap frequencies, although they use the same Newick field
- **<a id="gloss-use-11"></a>Transfer bootstrap expectation <sup>[11](#gloss-11)</sup>** (TBE) <a id="cite-57"></a>[Lemoine et al. 2018](https://doi.org/10.1038/s41586-018-0043-0) [[57](#ref-57)] measures how many taxa must move to recover a branch in each replicate tree, which keeps deep branches of large trees informative
- **Bayesian bootstrap** gives a single-mutation branch an expected support of about 90%, against about 63% for the Felsenstein bootstrap <a id="cite-10b"></a>[Lemoine and Gascuel 2024](https://doi.org/10.1093/molbev/msae238) [[10](#ref-10)]
- **<a id="gloss-use-12"></a>SPRTA <sup>[12](#gloss-12)</sup>** <a id="cite-11b"></a>[De Maio et al. 2025](https://doi.org/10.1038/s41586-025-09567-x) [[11](#ref-11)] scores each branch $b$ by the probability that its lower node descends directly from its upper node, against the alternative attachments of the subtree $S_b$:

$$
\mathrm{SPRTA}(b) = \frac{P(D \mid T)}{\sum_{i=1}^{I_b} P(D \mid T_i^b)}
$$

where $D$ is the alignment and $T_i^b$, for $i = 1, \dots, I_b$, are the topologies obtained by regrafting $S_b$ at alternative positions, with $T_1^b = T$. SPRTA evaluates the same SPR moves that the tree search evaluates, so it can run during the search "at negligible additional computational cost"; it also scores terminal branches (the placement of single samples), and is at least two orders of magnitude cheaper than bootstrap methods. Run after the search on a tree of 2,072,111 genomes, it took 7 h 27 min and 26.93 GB on one core

- **Parsimony placement uncertainty.** UShER reports the number of equally parsimonious placements of each sample <a id="cite-26c"></a>[Turakhia et al. 2021](https://doi.org/10.1038/s41588-021-00862-7) [[26](#ref-26)], and matUtils summarizes it per sample <a id="cite-28b"></a>[McBroome et al. 2021](https://doi.org/10.1093/molbev/msab264) [[28](#ref-28)]

TreeTime v1 reads input support values and writes them nowhere, because the runs reroot, collapse and resolve branches and the values must stay with their split ([kb/issues/M-io-branch-support-dropped-from-outputs.md](../issues/M-io-branch-support-dropped-from-outputs.md)). Support that TreeTime computes itself has the same requirement.

### Time-aware topology inference

- **Bayesian joint inference.** BEAST X <a id="cite-58"></a>[Baele et al. 2025](https://doi.org/10.1038/s41592-025-02751-x) [[58](#ref-58)] samples topologies, dates, and model parameters together and adds gradient-informed samplers for high-dimensional parameters. MCMC over timed trees is "limited to around a few thousand sequences unless tree space is restricted"; parsimony-informed topology operators reduce the time to convergence <a id="cite-59"></a>[Bouckaert et al. 2025](https://doi.org/10.1101/2025.06.18.660471) [[59](#ref-59)]. Online BEAST inserts new sequences into a converged analysis and resumes it, which shortens the time to an updated posterior <a id="cite-60"></a>[Gill et al. 2020](https://doi.org/10.1093/molbev/msaa047) [[60](#ref-60)]
- **Delphy** <a id="cite-31b"></a>[Varilly et al. 2026](https://doi.org/10.1038/s41586-026-11012-6) [[31](#ref-31)] reformulates Bayesian phylogenetics on an <a id="gloss-use-13"></a>explicit mutation-annotated tree <sup>[13](#gloss-13)</sup>: a timed tree whose branches carry the individual mutations with their times. Moves cost time in proportion to the mutations near the changed branches, and a mutation-directed SPR move proposes regrafting points from sequence differences. The paper reports 2 to 3 orders of magnitude speedups over the previous state of the art and a 100,000-sequence analysis within a day. Delphy has the MIT license since 2025-12-17 (release 1.4.1, 2026-07-15) and builds its own starting tree (see [Parsimony and placement on mutation-annotated trees](#parsimony-and-placement-on-mutation-annotated-trees))
- **Effect of a fixed topology.** <a id="cite-4b"></a>[Fourment et al. 2026](https://doi.org/10.1093/sysbio/syag069) [[4](#ref-4)] measured the bias of the two-step approach (see [Why the topology stage matters for dating](#why-the-topology-stage-matters-for-dating)). Posterior tree landscapes of large datasets are "diffuse yet rugged", and a small number of sequences such as recombinants and recurrent mutants causes most mixing problems <a id="cite-61"></a>[Gao et al. 2026](https://doi.org/10.1073/pnas.2510938123) [[61](#ref-61)]
- **Dating on a fixed topology.** LSD <a id="cite-62"></a>[To et al. 2016](https://doi.org/10.1093/sysbio/syv068) [[62](#ref-62)] fits node dates by least squares in linear time on a rooted tree. Chronumental <a id="cite-63"></a>[Sanderson 2021](https://doi.org/10.1101/2021.10.27.465994) [[63](#ref-63)] fits dates by stochastic gradient descent and handles trees with millions of nodes. Both take the topology as input and are comparators for TreeTime's dating, not tree builders
- **Gradients.** Branch-length gradients of the tree likelihood can be computed for all branches in linear time <a id="cite-64"></a>[Ji et al. 2020](https://doi.org/10.1093/molbev/msaa130) [[64](#ref-64)]. The same pre-order and post-order passes that TreeTime runs for marginal reconstruction provide the quantities these gradients need
- **TreeTime's polytomy resolution** groups children of a polytomy pairwise by the gain in temporal likelihood <a id="cite-1c"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[1](#ref-1)]. It is a local topology inference step that uses dates. This survey found no published ML tree search that scores topology moves with tip dates. A time-aware SPR search, scored by the sum of the sequence log-likelihood and TreeTime's temporal log-likelihood, would be a new method that needs its own validation

### Recombination and ancestral recombination graphs

A single tree assumes that all sites share one history. Methods that relax this:

- **RIPPLES** <a id="cite-65"></a>[Turakhia et al. 2022](https://doi.org/10.1038/s41586-022-05189-9) [[65](#ref-65)] splits candidate nodes of a MAT at one or two breakpoints and tests whether partial placements reduce the parsimony score
- **sc2ts** <a id="cite-66"></a>[Zhan et al. 2023](https://doi.org/10.1101/2023.06.08.544212) [[66](#ref-66)] (preprint) infers an <a id="gloss-use-14"></a>ancestral recombination graph <sup>[14](#gloss-14)</sup> (ARG) for 2.48 million SARS-CoV-2 genomes in daily batches with a copying model from the tsinfer framework <a id="cite-67"></a>[Kelleher et al. 2019](https://doi.org/10.1038/s41588-019-0483-y) [[67](#ref-67)]
- **Gubbins** <a id="cite-68"></a>[Croucher et al. 2015](https://doi.org/10.1093/nar/gku1196) [[68](#ref-68)] and **ClonalFrameML** <a id="cite-69"></a>[Didelot and Wilson 2015](https://doi.org/10.1371/journal.pcbi.1004041) [[69](#ref-69)] detect recombined regions in bacterial genomes and infer the clonal tree from the remaining sites

None of these tools appears in the inspected Nextstrain workflows. TreeTime v1 traversals accept only one parent per node (`fn Graph.one_parent_of()` [`packages/treetime-graph/src/graph.rs#L63-L74`](../../packages/treetime-graph/src/graph.rs#L63-L74)), and the v0 two-tree ARG features are not present in v1 ([kb/features/arg.md](../features/arg.md)). Recombination-aware inference changes the data model and the output contract, so it is a separate product decision from tree building.

### Machine learning

- **End-to-end builders.** Phyloformer <a id="cite-70"></a>[Nesterenko et al. 2025](https://doi.org/10.1093/molbev/msaf051) [[70](#ref-70)] predicts pairwise distances with a neural network and builds the tree with FastME. It was trained and evaluated on protein models, needs a GPU for its speed advantage, and its topological accuracy falls behind ML methods as the number of sequences grows. NeuralNJ <a id="cite-71"></a>[Zhang et al. 2025](https://doi.org/10.1093/molbev/msaf260) [[71](#ref-71)] learns the join steps of NJ, and IQ-NET <a id="cite-72"></a>[Yang et al. 2027](https://doi.org/10.1016/j.ympev.2026.108744) [[72](#ref-72)] infers quartet topologies with a network trained on empirical alignments
- **Search guidance.** A trained model can rank SPR candidates <a id="cite-73"></a>[Azouri et al. 2021](https://doi.org/10.1038/s41467-021-22073-8) [[73](#ref-73)], and reinforcement learning can drive the search on small datasets <a id="cite-74"></a>[Azouri et al. 2024](https://doi.org/10.1093/molbev/msae105) [[74](#ref-74)]
- **Production use.** RAxML-NG uses difficulty prediction to set its search effort <a id="cite-38b"></a>[Haag et al. 2022](https://doi.org/10.1093/molbev/msac254) [[38](#ref-38)], and IQ-TREE 3 ships neural-network models for model selection (`nn_models/` in the IQ-TREE 3 source)
- **Automatic differentiation.** General-purpose automatic differentiation computes phylogenetic gradients much more slowly than dedicated gradient code <a id="cite-75"></a>[Fourment et al. 2023](https://doi.org/10.1093/gbe/evad099) [[75](#ref-75)]

**Fit for TreeTime.** No machine-learning builder in this survey supports nucleotide outbreak alignments with thousands of sequences on a CPU. Difficulty prediction and move ranking are optional accelerators for a search that already has a defined objective.

## Software ecosystem and licenses

TreeTime is MIT-licensed. The Free Software Foundation's GPL FAQ states that a program that links a GPL library must be licensed under the GPL as a whole [[doc](https://www.gnu.org/licenses/gpl-faq.html#IfLibraryIsGPL)], and that "pipes, sockets and command-line arguments are communication mechanisms normally used between two separate programs" [[doc](https://www.gnu.org/licenses/gpl-faq.html#MereAggregation)]. The lists below group the tools by what this permits.

**GPL or AGPL: separate process, or algorithm reference for a reimplementation**

- IQ-TREE 3, v3.1.4 (2026-09-10), GPL-2.0 [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/LICENSE)]. A static library mode exposes `build_tree`, `fit_tree`, `modelfinder`, and `build_njtree` as C functions [[src](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/main/libiqtree_fun.h#L73-L110)]; linking it makes the combined program GPL
- CMAPLE, v2.0.0 (2026-05-01), GPL-2.0; MAPLE, v0.7.5 (2025-10-16), GPL-3.0
- RAxML-NG, 2.0.3 (2026-09-03), AGPL-3.0 [[src](https://github.com/amkozlov/raxml-ng/blob/d396351ee1263b704adb135edd6a4a84522dbd5b/LICENSE.txt)], and its kernel library coraxlib, AGPL-3.0. The AGPL also covers network use, which matters for the TreeTime web server
- FastTree 2.2.0 (2025-06-02), GPL (the repository license file is GPL-3.0); VeryFastTree v4.0.5 (2025-04-06), GPL-3.0
- DecentTree v1.0.0 (2023-11-05), GPL-2.0; RapidNJ, GPL-2.0; FastME, GPL-3.0

**Permissive: can be embedded, ported, or read as reference code**

- UShER and matOptimize, v0.6.6 (2025-07-11, pre-release), MIT [[src](https://github.com/yatisht/usher/blob/ac9c982d937c3bc43e2c16a0b73d30cf2b937118/LICENSE)]
- Nextclade, 3.24.0 (2026-09-30), MIT, written in Rust and compiled to WebAssembly for its web application [[src](https://github.com/nextstrain/nextclade/blob/529db53e8f5b16d8a5d1161abd8b60ef04109817/LICENSE)]
- Delphy, 1.4.1 (2026-07-15), MIT [[src](https://github.com/broadinstitute/delphy/blob/5d5989f3e4cac9ed5ae9bbcad62bcbc81afd8c2a/LICENSE)]
- BEAGLE, v4.0.1 (2023-10-13), MIT, a likelihood computation library without tree search [[doc](https://github.com/beagle-dev/beagle-lib)]
- Rust crates, all young or small: `phyne` (ML and parsimony search with SPR and NNI moves, unpublished on crates.io, v0.1.2) [[doc](https://github.com/acg-team/phyne)], `ninja-phylo` (exact large-scale NJ, BSD-3-Clause, 2.0.0-rc.4) [[doc](https://github.com/TravisWheelerLab/ninja)], `nj` (MIT) [[doc](https://github.com/holmrenser/nj.rs)], `phylo` (tree structures and likelihood without search, MIT) [[doc](https://github.com/sriram98v/phylo-rs)]

**Consequences**

- No permissively licensed production ML search engine exists. A TreeTime-internal builder is a reimplementation from published algorithms, with the GPL tools as external test oracles
- The low-divergence references (UShER, matOptimize, Nextclade, Delphy) are MIT-licensed, so their code can be read and ported
- The TreeTime web and desktop applications need every component in-process or bundled. A GPL builder shipped next to TreeTime stays a separate program with its own license obligations

## TreeTime v1 building blocks and gaps

### Components that a builder can reuse

- **Graph edits.** `fn Graph.add_node()` [`packages/treetime-graph/src/graph_ops.rs#L12`](../../packages/treetime-graph/src/graph_ops.rs#L12), `fn Graph.add_edge()` [`packages/treetime-graph/src/graph_ops.rs#L47`](../../packages/treetime-graph/src/graph_ops.rs#L47), `fn Graph.reparent_edge()` [`packages/treetime-graph/src/graph_ops.rs#L93`](../../packages/treetime-graph/src/graph_ops.rs#L93), and `fn Graph.remove_edge()` [`packages/treetime-graph/src/graph_ops.rs#L143`](../../packages/treetime-graph/src/graph_ops.rs#L143) are the primitives of insertion and SPR moves; `fn MarginalReconstruction.reconcile_topology()` [`packages/treetime/src/partition/marginal/reconstruction.rs#L286`](../../packages/treetime/src/partition/marginal/reconstruction.rs#L286) adapts partition data to a changed topology
- **Fitch compression.** `fn create_fitch_partition()` [`packages/treetime/src/partition/fitch/passes.rs#L25`](../../packages/treetime/src/partition/fitch/passes.rs#L25) computes the root sequence and the substitutions on each edge, the content of a MAT. `fn infer_gtr_fitch()` [`packages/treetime/src/partition/fitch/gtr_inference.rs#L12`](../../packages/treetime/src/partition/fitch/gtr_inference.rs#L12) estimates an initial GTR model from the Fitch mutation counts
- **Sparse marginal likelihood.** Each edge stores the message from the child subtree (`struct SparseEdgeBackward` with `msg_to_parent` and `msg_from_child`, [`packages/treetime/src/partition/storage/sparse.rs#L141`](../../packages/treetime/src/partition/storage/sparse.rs#L141)) and the message from the rest of the tree (`struct SparseEdgeForward` with `msg_to_child` and `msg_from_parent`, [`packages/treetime/src/partition/storage/sparse.rs#L147`](../../packages/treetime/src/partition/storage/sparse.rs#L147)). These are the <a id="gloss-use-15"></a>inside and outside messages <sup>[15](#gloss-15)</sup> that placement and rearrangement scores need (see [Placement and rearrangement scores](#placement-and-rearrangement-scores-from-existing-messages))
- **Dense marginal likelihood.** `struct PartitionMarginalDense` [`packages/treetime/src/partition/marginal/dense/partition.rs#L32`](../../packages/treetime/src/partition/marginal/dense/partition.rs#L32) gives an independent implementation for comparison with the sparse one
- **Branch-length optimization.** `fn run_optimize_mixed()` [`packages/treetime/src/optimize/dispatch.rs#L18`](../../packages/treetime/src/optimize/dispatch.rs#L18) optimizes all edges from their messages with Newton or Brent methods, including an indel term
- **Parsimony-guided topology moves.** The `optimize` loop collapses internal edges whose optimal length is zero (`fn find_zero_optimal_internal_edges()` [`packages/treetime/src/optimize/run_loop.rs#L236`](../../packages/treetime/src/optimize/run_loop.rs#L236)) and resolves polytomies (`fn resolve_polytomies()` [`packages/treetime/src/optimize/topology/resolve_polytomy.rs#L18`](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs#L18)): it groups siblings that share substitutions (`fn merge_single_polytomy()` [`packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L39`](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs#L39)) and moves a child that reverts a substitution of its parent edge (`fn hoist_reverting_child()` [`packages/treetime/src/optimize/topology/hoist_reversions.rs#L132`](../../packages/treetime/src/optimize/topology/hoist_reversions.rs#L132)). These moves cover the operations of the mpox `fix_tree` script. [kb/decisions/optimize-polytomy-reversion-resolution.md](../decisions/optimize-polytomy-reversion-resolution.md) records them and names subtree parsimony re-optimization as the exact method, outside the scope of that decision
- **Rerooting.** `fn reroot_min_dev()` [`packages/treetime/src/reroot/orchestrate.rs#L20`](../../packages/treetime/src/reroot/orchestrate.rs#L20) roots by root-to-tip regression against dates, the same principle as Delphy's initializer, and `fn reroot_sparse()` [`packages/treetime/src/partition/marginal/sparse/reroot.rs#L13`](../../packages/treetime/src/partition/marginal/sparse/reroot.rs#L13) reuses the sparse messages after a reroot
- **Time-aware polytomy resolution.** `fn resolve_polytomies()` [`packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L16`](../../packages/treetime/src/timetree/optimization/polytomy/resolve.rs#L16) in the timetree pipeline
- **MAT input and output.** `packages/util-usher-mat` reads and writes UShER mutation-annotated trees ([kb/decisions/io-usher-mat-gaps-as-missing-data.md](../decisions/io-usher-mat-gaps-as-missing-data.md))

### Gaps

- **No starting topology.** Every command requires `--tree` (see [TreeTime v0 and v1](#treetime-v0-and-v1)); [kb/algo/unimplemented.md](../algo/unimplemented.md) lists tree inference as unimplemented
- **No global rearrangement search.** The topology moves are local to polytomies; there is no NNI or SPR search
- **No incremental update.** After a topology change, the `optimize` loop rebuilds the graph and recomputes all messages (`fn prune_and_merge_in_loop()` [`packages/treetime/src/optimize/run_loop.rs#L262`](../../packages/treetime/src/optimize/run_loop.rs#L262)). A search that evaluates many moves needs either local message updates or a cheaper surrogate score
- **Parsimony at multifurcations.** The Fitch backward pass uses the plurality recurrence at nodes with three or more children, so it keeps only minimum-cost states ([kb/decisions/ancestral-fitch-plurality-on-multifurcations.md](../decisions/ancestral-fitch-plurality-on-multifurcations.md)). A parsimony placement or SPR score can use it directly
- **Dense-sparse agreement.** [kb/issues/M-ancestral-dense-sparse-divergence.md](../issues/M-ancestral-dense-sparse-divergence.md) reports unexplained differences. A search that ranks topologies by likelihood can turn a small score difference into a different topology
- **No branch support,** and input support values are not written ([kb/issues/M-io-branch-support-dropped-from-outputs.md](../issues/M-io-branch-support-dropped-from-outputs.md))
- **No pairwise distances.** The only distance code is a Jukes-Cantor correction for merged branch lengths (`fn jukes_cantor_distance()` [`packages/treetime/src/gtr/jc_distance.rs#L5`](../../packages/treetime/src/gtr/jc_distance.rs#L5))

## Applicability to TreeTime

### Placement and rearrangement scores from existing messages

For an edge $e$ from parent $p$ to child $c$ with length $t_e$, let $f_{e,s}$ be the outside vector at $p$ for site $s$ (the message `msg_to_child`, which includes the root prior) and $g_{e,s}$ the inside vector at $c$ (the message `msg_to_parent`). Attaching a query sequence with observation vectors $h_s$ at a point that divides $e$ into lengths $t_1 + t_2 = t_e$, with a pendant branch of length $t_q$, gives

$$
\ell_e(t_1, t_q) = \sum_{s=1}^{L} \log \sum_{x} \big[P(t_1)^{\top} f_{e,s}\big]_x \, \big[P(t_2)\, g_{e,s}\big]_x \, \big[P(t_q)\, h_s\big]_x
$$

where $x$ runs over the states of the attachment point and $P(t) = e^{tQ}$. This is the three-taxon evaluation of pplacer (see [Placement on a fixed backbone tree](#placement-on-a-fixed-backbone-tree)). Consequences for TreeTime:

- After one full marginal pass, every candidate edge can be scored from stored messages, without recomputation
- In the sparse representation, the sum over sites splits into positions that are variable in $f_e$, in $g_e$, or in the query, which need explicit vectors, and the remaining positions, which group by state with a count. The cost per candidate edge grows with the number of variable positions near the edge and with the alphabet size, not with $L$
- An SPR move of a subtree $S$ replaces $h_s$ with the inside vector of $S$. The outside messages of the remaining tree also change when $S$ is removed, so a score from the current messages is an approximation; MAPLE and SPRTA accept candidates by such scores and then optimize the branch lengths of the survivors
- The parsimony counterpart replaces vectors by Fitch state sets and the score by the number of added substitutions, as in UShER and matOptimize

The code has the messages, but it has no function that scores a candidate attachment from them and no partial update after a change.

### Options

The choices below are independent unless a coupling is stated.

- **Entry point without `--tree`.** Keep the error; call an external builder as a subprocess, as v0 does; or run a native builder
- **Initial topology.** NJ or BIONJ on corrected distances; parsimony placement on the Fitch mutation representation, as in UShER, Nextclade, and Delphy's initializer; or likelihood placement on the sparse messages, as in MAPLE
- **Topology improvement.** None; parsimony SPR, as in matOptimize; sparse-likelihood SPR, as in MAPLE; or time-aware SPR. The improvement objective couples with the data structure of the initial topology: parsimony SPR works on the Fitch representation, likelihood SPR on the sparse messages
- **Branch support.** None; counts of equally parsimonious placements; SPRTA-style support, which needs likelihood SPR scores; or bootstrap variants
- **Update mode.** Build from scratch in each run, or place new sequences onto a saved tree and refine locally. Incremental update requires placement
- **Divergent data.** Restrict the native builder to data that passes a divergence test, as IQ-TREE's `--pathogen` does, and require an input tree otherwise; or add a distance tree with likelihood NNI refinement for divergent data

### Assessment against the evidence

- **Parsimony placement with parsimony SPR** has the strongest evidence for the low-divergence datasets, uses the representation TreeTime already builds, and has MIT-licensed reference implementations. Its output is a multifurcating tree with mutations on edges, which feeds the existing `optimize` loop (ML branch lengths, shared-mutation merge, reversion hoist) and the timetree polytomy resolution directly. It depends on a correct multifurcation recurrence and on deterministic tie handling. Placement order matters at scale; TreeTime has sampling dates, which allow the temporal order that the Viridian tree used
- **NJ or BIONJ** works on any divergence and gives a deterministic binary tree. At low divergence it resolves short branches arbitrarily, and its distance matrix needed about 50 GB for 64,000 sequences in DecentTree unless the sparse or dynamic NJ variants are used
- **Sparse-likelihood SPR** follows MAPLE and CMAPLE, which have the highest reported accuracy at low divergence. It needs a candidate-scoring function on the stored messages, local message updates, and resolved dense-sparse agreement. CMAPLE's tree code has about 15,000 lines of C++, a measure of the code size such a feature adds
- **Time-aware SPR** targets the bias described in [Why the topology stage matters for dating](#why-the-topology-stage-matters-for-dating), and TreeTime holds both the sequence and the temporal likelihood. No published method exists to compare against, so it is a research project with its own validation
- **External subprocess fallback** restores v0 behavior and keeps IQ-TREE's accuracy on divergent data, but it does not meet the goal of operation without external tools and is unavailable in the web and desktop applications unless a builder is bundled
- **Machine-learning builders** are not applicable now (see [Machine learning](#machine-learning))

### Validation approach

A native builder changes scientific outputs, so each part needs evidence of the kinds below before adoption:

- **Independent oracles.** Compare log-likelihoods under the same model and the same branch-length optimization with IQ-TREE, RAxML-NG, and CMAPLE run outside the TreeTime test suite, and store the results as reference fixtures. v0 cannot serve as an oracle, because it has no builder of its own
- **Simulated data with known trees,** at several divergence levels, scored by Robinson-Foulds and quartet distances after collapsing zero-length branches in both trees, because binary resolution of true polytomies is arbitrary
- **Representation checks.** Placement and SPR scores from stored messages must equal the log-likelihood of a full pass on the changed tree, separately for dense and sparse partitions
- **Downstream effect.** Compare clock rate, root date, and node-date intervals of the timetree pipeline with an input tree from IQ-TREE and with the native tree, on each low-divergence family in `data/`
- **Runtime and memory scaling** on the largest datasets of each family

## Open decisions

### Settled

- **Embedding a GPL or AGPL engine in TreeTime is excluded.** Linking makes the whole program GPL, by the GPL FAQ (see [Software ecosystem and licenses](#software-ecosystem-and-licenses))
- **Machine-learning builders do not support TreeTime's data now.** The published builders target protein data, small trees, or GPUs (see [Machine learning](#machine-learning))
- **Low-divergence data favors parsimony placement and sparse likelihood methods.** UShER with matOptimize equals or exceeds ML builders on SARS-CoV-2, MAPLE and CMAPLE exceed both, and the advantage decreases with divergence (see [Parsimony and placement on mutation-annotated trees](#parsimony-and-placement-on-mutation-annotated-trees) and [Pandemic-scale likelihood](#pandemic-scale-likelihood-maple-and-cmaple))
- **The v1 sparse likelihood keeps the constant-site contribution.** The fixed-state counts enter the log-likelihood, so an ascertainment correction is needed only for inputs without invariant sites (code evidence in [Pandemic-scale likelihood](#pandemic-scale-likelihood-maple-and-cmaple))
- **v0 has no internal tree builder.** It calls external programs; the KB issue that describes neighbor joining in v0 is out of date (see [TreeTime v0 and v1](#treetime-v0-and-v1))
- **Placement cannot build the first tree.** It needs a backbone, so a native builder needs an initial-topology method in addition to placement

### Open

Seven decisions remain. A native builder diverges from v0, which delegates to external programs, so each adopted option needs an approved entry in `kb/decisions/`.

#### Behavior without an input tree

Today every command stops with an error when `--tree` is missing. What should TreeTime do?

- **Keep the error.** Example: the user runs `augur tree` first, as today
- **Subprocess fallback.** Example: TreeTime calls `iqtree3` when it is installed, as v0 calls IQ-TREE, FastTree, or RAxML
- **[recommended] Native builder behind a divergence test.** Example: `treetime timetree` without `--tree` builds a tree for an mpox alignment and stops with an error that names an external builder for a dengue alignment

#### Initial topology

Which method builds the first tree?

- **NJ or BIONJ on corrected distances.** Example: a BIONJ tree from TN93 distances with pairwise deletion of gaps
- **[recommended] Parsimony placement on the Fitch mutation representation.** Example: insert samples in date order at the edge with the fewest added substitutions, as UShER and Delphy's initializer do
- **Likelihood placement on sparse messages.** Example: MAPLE-style placement that scores each candidate edge with the formula in [Placement and rearrangement scores](#placement-and-rearrangement-scores-from-existing-messages)

#### Topology improvement (multi-select)

Which rearrangement search improves the first tree? Options can be combined in sequence.

- **[recommended] Parsimony SPR, then the existing `optimize` loop.** Example: matOptimize-style SPR with a growing radius, then ML branch lengths and polytomy moves
- **Sparse-likelihood SPR.** Example: MAPLE-style SPR with local message updates. Requires resolved dense-sparse agreement
- **Time-aware SPR.** Example: SPR scored by the sum of the sequence and temporal log-likelihoods. Research without a published reference

#### Divergent data

What happens for alignments that fail the divergence test?

- **[recommended] Require an input tree.** Example: the error message for a Lassa alignment names IQ-TREE and the divergence values
- **Distance tree with likelihood NNI.** Example: BIONJ start, then NNI moves scored on dense messages

#### Branch support (multi-select)

Which support values should a native builder report? Options can be combined.

- **[recommended] Counts of equally parsimonious placements.** Example: each sample gets the number of edges with the same parsimony cost, as in UShER
- **SPRTA-style support.** Example: the probability of each branch against alternative regrafts. Requires likelihood SPR
- **Bayesian or transfer bootstrap.** Example: 100 replicate searches summarized as Bayesian bootstrap support

#### Update mode

Should TreeTime accept a saved tree plus new sequences?

- **Build from scratch only.** Example: every run builds the full tree
- **[recommended] Also place new sequences onto a saved tree.** Example: a MAT from the previous run plus 500 new sequences, placed and refined locally

#### Recombination

Should the builder handle recombination?

- **[recommended] Outside the builder's scope.** Example: the documentation states the single-tree assumption, and masked sites stay a pipeline step
- **Accept a recombination mask.** Example: a BED file of regions to exclude, as Gubbins produces
- **ARG inference.** Example: sc2ts-style copying model. A separate product objective

#### Recommended combination

A native builder behind a divergence test, with parsimony placement in date order as the initial topology, parsimony SPR followed by the existing `optimize` loop and timetree polytomy resolution, equally parsimonious placement counts as support, incremental placement onto saved trees, an input tree required for divergent data, and recombination outside the scope. This combination reuses the Fitch representation, the topology moves, and the MAT input and output of TreeTime v1, and the evidence supports it for the low-divergence datasets. Sparse-likelihood SPR and time-aware SPR are the next steps once the scoring function and local updates exist. The Fitch multifurcation recurrence and the dense-sparse divergence need fixes before either score can rank topologies.

## Cross-topic themes

- **Difference-based representations set the scale.** UShER's MAT, MAPLE's genome lists, Delphy's explicit mutation-annotated tree, and TreeTime's sparse partitions all store sequences as differences from a reference or an ancestor. Their cost grows with the number of mutations, not with the genome length
- **Placement with local refinement replaces global search at scale.** UShER with matOptimize, MAPLE, CMAPLE, Nextclade, and Delphy's initializer all insert samples one at a time and then repair the tree locally
- **Polytomies are the correct output at low divergence.** A binary tree forces an arbitrary order on unresolved branchings. Methods that keep zero-length branches collapsed (UShER, multifurcating NJ, IQ-TREE `--polytomy`) give TreeTime's polytomy resolution its input directly
- **Support moves from split frequencies to placement probabilities.** SPRTA and counts of equally parsimonious placements describe where a lineage attaches, which matches the questions of genomic epidemiology
- **Dates enter topology inference rarely.** Only Bayesian tools (BEAST, Delphy) and local heuristics (TreeTime's polytomy resolution, the temporal insertion order of the Viridian tree) use sampling dates while they choose the topology

## Emerging trends

- IQ-TREE 3 (2025-2026) selects CMAPLE automatically for low-divergence data and adds SPRTA
- RAxML-NG 2.0 (2026) makes the difficulty-adaptive search the default and adds early stopping and machine-learning branch support prediction
- The MAPLE line adds rate and sequencing-error models and local references (2026), and multithreaded search in CMAPLE 2 (2026 preprint)
- Delphy (2026) brings Bayesian joint inference of topology and dates to 100,000 sequences, with an MIT license and a browser application
- Trees of several million samples, built by placement and parsimony optimization, show that insertion order and systematic sequencing errors control the deep structure of the tree

## Controversies and conflicting evidence

- **Parsimony against ML at low divergence.** The preprint of the online-parsimony study reported "slightly better trees" for parsimony <a id="cite-76"></a>[Thornlow et al. 2021](https://doi.org/10.1101/2021.12.02.471004) [[76](#ref-76)]; the published version reports "equivalent trees" <a id="cite-9c"></a>[Kramer et al. 2023](https://doi.org/10.1093/sysbio/syad031) [[9](#ref-9)]. In the MAPLE benchmarks, matOptimize was less accurate than ML methods on simulated data and more accurate on real data. The conclusion depends on whether simulated or real data is the reference
- **Rate model.** Model selection often chooses the FreeRate model, while <a id="cite-8c"></a>[Morel et al. 2021](https://doi.org/10.1093/molbev/msaa314) [[8](#ref-8)] recommend the discrete gamma model for SARS-CoV-2 because FreeRate likelihoods were numerically unstable in the search
- **Bootstrap calibration.** The Felsenstein bootstrap underestimates the support of short correct branches, and the ultrafast bootstrap needed a correction for overestimated support on polytomies. No support method is calibrated for both cases
- **Reference for the fixed-topology bias.** <a id="cite-4c"></a>[Fourment et al. 2026](https://doi.org/10.1093/sysbio/syag069) [[4](#ref-4)] measure bias against an unconstrained Bayesian analysis, which is itself model-based, and use TreeTime only for rooting
- **Stated licenses.** The RAxML-NG 2.0 preprint states that the code is available "under GNU GPL"; the repository license file is the GNU Affero GPL version 3
- **matOptimize stopping rule.** The paper states 0.5% and the code default is 0.05%

## Gaps and open questions

- **Evidence base.** The parsimony-against-ML comparisons are mostly on SARS-CoV-2. Other low-divergence pathogens in `data/` (mpox, Ebola, RSV) have no comparable published benchmark
- **Time-aware search.** This survey found no ML tree search that uses tip dates during the topology search, and no measurement of how a polytomy-preserving builder changes TreeTime's dating
- **Partly read sources.** MAPLE's supplementary methods (the exact genome-list rules) were not read; the rules are visible in the CMAPLE source. The Delphy statements come from the abstract and the source code. The PhyML, EPA-ng, APPLES, SCAMPP, Gubbins, ClonalFrameML, BEAST X, and TargetedBeast statements come from abstracts
- **Dataset divergence.** The H3N2 HA set fails the CMAPLE test by a small margin (mean 0.0222 against 0.02). The full-tree influenza builds force CMAPLE regardless of the test

## Glossary

1. <a id="gloss-1"></a> **Topology.** The branching pattern of a tree: which samples and ancestors are grouped together, independent of branch lengths. [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Polytomy.** An internal node with more than two children. It represents either simultaneous divergence or a branching order that the data cannot resolve. [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **Ascertainment bias.** The distortion of likelihood estimates when the alignment contains only sites selected for variation; an ascertainment correction conditions the likelihood on the selection rule. [↩](#gloss-use-3)
4. <a id="gloss-4"></a> **Balanced minimum evolution (BME).** A distance criterion that scores a tree by its total branch length, estimated with a weighting that gives the two subtrees at each internal node equal weight; NJ optimizes it greedily ([Gascuel and Steel 2006](https://doi.org/10.1093/molbev/msl072) [[13](#ref-13)]). [↩](#gloss-use-4)
5. <a id="gloss-5"></a> **Mutation-annotated tree (MAT).** A tree that stores a reference or root sequence and, on each branch, the list of mutations along that branch. UShER's protobuf format is the common file format. [↩](#gloss-use-5)
6. <a id="gloss-6"></a> **Subtree pruning and regrafting (SPR).** A topology move that detaches a subtree and attaches it to another edge. The SPR radius limits how far, in edges, the subtree can move. [↩](#gloss-use-6)
7. <a id="gloss-7"></a> **Long-branch attraction.** The tendency of parsimony, and of misspecified models, to group long branches together because parallel changes on them look like shared ancestry ([Felsenstein 1978](https://doi.org/10.1093/sysbio/27.4.401) [[32](#ref-32)]). [↩](#gloss-use-7)
8. <a id="gloss-8"></a> **Nearest-neighbor interchange (NNI).** A topology move that swaps two subtrees across one internal edge; each internal edge has two NNI alternatives. [↩](#gloss-use-8)
9. <a id="gloss-9"></a> **Genome list.** MAPLE's representation of a partial likelihood: run-length entries that cover runs of reference positions with one entry and store explicit vectors only for differing or uncertain positions. [↩](#gloss-use-9)
10. <a id="gloss-10"></a> **Phylogenetic placement.** Attaching new sequences to the edges of an existing tree, the backbone, without changing the backbone. [↩](#gloss-use-10)
11. <a id="gloss-11"></a> **Transfer bootstrap expectation (TBE).** A bootstrap support that counts, for each replicate tree, the minimum number of taxa that must move to recover a branch, divides it by one less than the number of taxa on the smaller side of the branch, and averages one minus this value over the replicates. [↩](#gloss-use-11)
12. <a id="gloss-12"></a> **SPRTA.** Subtree pruning and regrafting-based tree assessment: a branch support equal to the probability that the lower node of a branch descends directly from its upper node, relative to alternative regrafts. [↩](#gloss-use-12)
13. <a id="gloss-13"></a> **Explicit mutation-annotated tree.** Delphy's timed tree in which each branch carries the individual mutations with their times, so the likelihood factorizes into terms near each mutation. [↩](#gloss-use-13)
14. <a id="gloss-14"></a> **Ancestral recombination graph (ARG).** A graph of ancestry in which a node can have two parents, so that different genome segments follow different trees. [↩](#gloss-use-14)
15. <a id="gloss-15"></a> **Inside and outside messages.** For an edge, the inside message is the partial likelihood of the data in the subtree below the edge, and the outside message is the partial likelihood of all other data, including the root prior. Their product over an edge gives the full likelihood. [↩](#gloss-use-15)

## References

1. <a id="ref-1"></a> Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-likelihood phylodynamic analysis." _Virus Evolution_ 4 (1): vex042. https://doi.org/10.1093/ve/vex042 [↩¹](#cite-1a) [↩²](#cite-1b) [↩³](#cite-1c)
2. <a id="ref-2"></a> Hadfield, James, Colin Megill, Sidney M. Bell, et al. 2018. "Nextstrain: Real-time tracking of pathogen evolution." _Bioinformatics_ 34 (23): 4121-4123. https://doi.org/10.1093/bioinformatics/bty407 [↩](#cite-2)
3. <a id="ref-3"></a> Huddleston, John, James Hadfield, Thomas Sibley, et al. 2021. "Augur: A bioinformatics toolkit for phylogenetic analyses of human pathogens." _Journal of Open Source Software_ 6 (57): 2906. https://doi.org/10.21105/joss.02906 [↩](#cite-3)
4. <a id="ref-4"></a> Fourment, Mathieu, Jiansi Gao, Marc A. Suchard, and Frederick A. Matsen IV. 2026. "Assessing the validity of the fixed tree topology assumption in phylodynamic inference." _Systematic Biology_, syag069. https://doi.org/10.1093/sysbio/syag069 [↩¹](#cite-4a) [↩²](#cite-4b) [↩³](#cite-4c)
5. <a id="ref-5"></a> Wong, Thomas K. F., Nhan Ly-Trong, Huaiyan Ren, et al. 2026. "IQ-TREE 3: Phylogenomic inference software using complex evolutionary models." _Molecular Biology and Evolution_ 43 (5): msag117. https://doi.org/10.1093/molbev/msag117 [↩¹](#cite-5a) [↩²](#cite-5b)
6. <a id="ref-6"></a> Ly-Trong, Nhan, Chris Bielow, Nicola De Maio, and Bui Quang Minh. 2024. "CMAPLE: Efficient phylogenetic inference in the pandemic era." _Molecular Biology and Evolution_ 41 (7): msae134. https://doi.org/10.1093/molbev/msae134 [↩¹](#cite-6a) [↩²](#cite-6b)
7. <a id="ref-7"></a> Wang, Weiwen, James Barbetti, Thomas Wong, et al. 2023. "DecentTree: Scalable neighbour-joining for the genomic era." _Bioinformatics_ 39 (9): btad536. https://doi.org/10.1093/bioinformatics/btad536 [↩¹](#cite-7a) [↩²](#cite-7b)
8. <a id="ref-8"></a> Morel, Benoit, Pierre Barbera, Lucas Czech, et al. 2021. "Phylogenetic analysis of SARS-CoV-2 data is difficult." _Molecular Biology and Evolution_ 38 (5): 1777-1791. https://doi.org/10.1093/molbev/msaa314 [↩¹](#cite-8a) [↩²](#cite-8b) [↩³](#cite-8c)
9. <a id="ref-9"></a> Kramer, Alexander M., Bryan Thornlow, Cheng Ye, et al. 2023. "Online phylogenetics with matOptimize produces equivalent trees and is dramatically more efficient for large SARS-CoV-2 phylogenies than de novo and maximum-likelihood implementations." _Systematic Biology_ 72 (5): 1039-1051. https://doi.org/10.1093/sysbio/syad031 [↩¹](#cite-9a) [↩²](#cite-9b) [↩³](#cite-9c)
10. <a id="ref-10"></a> Lemoine, Frédéric, and Olivier Gascuel. 2024. "The Bayesian phylogenetic bootstrap and its application to short trees and branches." _Molecular Biology and Evolution_ 41 (11): msae238. https://doi.org/10.1093/molbev/msae238 [↩¹](#cite-10a) [↩²](#cite-10b)
11. <a id="ref-11"></a> De Maio, Nicola, Nhan Ly-Trong, Samuel Martin, Bui Quang Minh, and Nick Goldman. 2025. "Assessing phylogenetic confidence at pandemic scales." _Nature_ 647 (8089): 472-478. https://doi.org/10.1038/s41586-025-09567-x [↩¹](#cite-11a) [↩²](#cite-11b)
12. <a id="ref-12"></a> Saitou, Naruya, and M. Nei. 1987. "The neighbor-joining method: A new method for reconstructing phylogenetic trees." _Molecular Biology and Evolution_ 4 (4): 406-425. https://doi.org/10.1093/oxfordjournals.molbev.a040454 [↩](#cite-12)
13. <a id="ref-13"></a> Gascuel, Olivier, and M. Steel. 2006. "Neighbor-joining revealed." _Molecular Biology and Evolution_ 23 (11): 1997-2000. https://doi.org/10.1093/molbev/msl072 [↩](#cite-13)
14. <a id="ref-14"></a> Gascuel, Olivier. 1997. "BIONJ: An improved version of the NJ algorithm based on a simple model of sequence data." _Molecular Biology and Evolution_ 14 (7): 685-695. https://doi.org/10.1093/oxfordjournals.molbev.a025808 [↩](#cite-14)
15. <a id="ref-15"></a> Desper, Richard, and Olivier Gascuel. 2002. "Fast and accurate phylogeny reconstruction algorithms based on the minimum-evolution principle." _Journal of Computational Biology_ 9 (5): 687-705. https://doi.org/10.1089/106652702761034136 [↩](#cite-15)
16. <a id="ref-16"></a> Lefort, Vincent, Richard Desper, and Olivier Gascuel. 2015. "FastME 2.0: A comprehensive, accurate, and fast distance-based phylogeny inference program." _Molecular Biology and Evolution_ 32 (10): 2798-2800. https://doi.org/10.1093/molbev/msv150 [↩](#cite-16)
17. <a id="ref-17"></a> Simonsen, Martin, Thomas Mailund, and Christian N. S. Pedersen. 2008. "Rapid neighbour-joining." In _Algorithms in Bioinformatics (WABI 2008)_, Lecture Notes in Computer Science, 113-122. Springer. https://doi.org/10.1007/978-3-540-87361-7_10 [↩](#cite-17)
18. <a id="ref-18"></a> Clausen, Philip T. L. C. 2023. "Scaling neighbor joining to one million taxa with dynamic and heuristic neighbor joining." _Bioinformatics_ 39 (1): btac774. https://doi.org/10.1093/bioinformatics/btac774 [↩](#cite-18)
19. <a id="ref-19"></a> Kurt, Semih, Alexandre Bouchard-Côté, and Jens Lagergren. 2024. "Sparse neighbor joining: Rapid phylogenetic inference using a sparse distance matrix." _Bioinformatics_ 40 (12): btae701. https://doi.org/10.1093/bioinformatics/btae701 [↩](#cite-19)
20. <a id="ref-20"></a> Price, Morgan N., Paramvir S. Dehal, and Adam P. Arkin. 2009. "FastTree: Computing large minimum evolution trees with profiles instead of a distance matrix." _Molecular Biology and Evolution_ 26 (7): 1641-1650. https://doi.org/10.1093/molbev/msp077 [↩](#cite-20)
21. <a id="ref-21"></a> Fernández, Alberto, Natàlia Segura-Alabart, and Francesc Serratosa. 2023. "The MultiFurcating neighbor-joining algorithm for reconstructing polytomic phylogenetic trees." _Journal of Molecular Evolution_ 91 (6): 773-779. https://doi.org/10.1007/s00239-023-10134-z [↩](#cite-21)
22. <a id="ref-22"></a> Tamura, Koichiro, and M. Nei. 1993. "Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees." _Molecular Biology and Evolution_ 10 (3): 512-526. https://doi.org/10.1093/oxfordjournals.molbev.a040023 [↩](#cite-22)
23. <a id="ref-23"></a> Fitch, W. M. 1971. "Toward defining the course of evolution: Minimum change for a specific tree topology." _Systematic Biology_ 20 (4): 406-416. https://doi.org/10.1093/sysbio/20.4.406 [↩](#cite-23)
24. <a id="ref-24"></a> Sankoff, David. 1975. "Minimal mutation trees of sequences." _SIAM Journal on Applied Mathematics_ 28 (1): 35-42. https://doi.org/10.1137/0128004 [↩](#cite-24)
25. <a id="ref-25"></a> Hartigan, John A. 1973. "Minimum mutation fits to a given tree." _Biometrics_ 29 (1): 53-65. https://doi.org/10.2307/2529676 [↩](#cite-25)
26. <a id="ref-26"></a> Turakhia, Yatish, Bryan Thornlow, Angie S. Hinrichs, et al. 2021. "Ultrafast sample placement on existing trees (UShER) enables real-time phylogenetics for the SARS-CoV-2 pandemic." _Nature Genetics_ 53 (6): 809-816. https://doi.org/10.1038/s41588-021-00862-7 [↩¹](#cite-26a) [↩²](#cite-26b) [↩³](#cite-26c)
27. <a id="ref-27"></a> Ye, Cheng, Bryan Thornlow, Angie Hinrichs, et al. 2022. "matOptimize: A parallel tree optimization method enables online phylogenetics for SARS-CoV-2." _Bioinformatics_ 38 (15): 3734-3740. https://doi.org/10.1093/bioinformatics/btac401 [↩](#cite-27)
28. <a id="ref-28"></a> McBroome, Jakob, Bryan Thornlow, Angie S. Hinrichs, et al. 2021. "A daily-updated database and tools for comprehensive SARS-CoV-2 mutation-annotated trees." _Molecular Biology and Evolution_ 38 (12): 5819-5824. https://doi.org/10.1093/molbev/msab264 [↩¹](#cite-28a) [↩²](#cite-28b)
29. <a id="ref-29"></a> Hunt, Martin, Angie S. Hinrichs, Daniel Anderson, et al. 2026. "Addressing pandemic-wide systematic errors in the SARS-CoV-2 phylogeny." _Nature Methods_ 23 (3): 653-662. https://doi.org/10.1038/s41592-025-02947-1 [↩](#cite-29)
30. <a id="ref-30"></a> Aksamentov, Ivan, Cornelius Roemer, Emma Hodcroft, and Richard Neher. 2021. "Nextclade: Clade assignment, mutation calling and quality control for viral genomes." _Journal of Open Source Software_ 6 (67): 3773. https://doi.org/10.21105/joss.03773 [↩](#cite-30)
31. <a id="ref-31"></a> Varilly, Patrick, Mark Schifferli, Katherine Yang, et al. 2026. "Scalable near-real-time Bayesian phylogenetics for outbreaks with Delphy." _Nature_. https://doi.org/10.1038/s41586-026-11012-6 [↩¹](#cite-31a) [↩²](#cite-31b)
32. <a id="ref-32"></a> Felsenstein, Joseph. 1978. "Cases in which parsimony or compatibility methods will be positively misleading." _Systematic Biology_ 27 (4): 401-410. https://doi.org/10.1093/sysbio/27.4.401 [↩](#cite-32)
33. <a id="ref-33"></a> De Maio, Nicola, Prabhav Kalaghatgi, Yatish Turakhia, Russell Corbett-Detig, Bui Quang Minh, and Nick Goldman. 2023. "Maximum likelihood pandemic-scale phylogenetics." _Nature Genetics_ 55 (5): 746-752. https://doi.org/10.1038/s41588-023-01368-0 [↩¹](#cite-33a) [↩²](#cite-33b)
34. <a id="ref-34"></a> Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh. 2015. "IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies." _Molecular Biology and Evolution_ 32 (1): 268-274. https://doi.org/10.1093/molbev/msu300 [↩](#cite-34)
35. <a id="ref-35"></a> Minh, Bui Quang, Heiko A. Schmidt, Olga Chernomor, et al. 2020. "IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era." _Molecular Biology and Evolution_ 37 (5): 1530-1534. https://doi.org/10.1093/molbev/msaa015 [↩](#cite-35)
36. <a id="ref-36"></a> Kozlov, Alexey M., Diego Darriba, Tomáš Flouri, Benoit Morel, and Alexandros Stamatakis. 2019. "RAxML-NG: A fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference." _Bioinformatics_ 35 (21): 4453-4455. https://doi.org/10.1093/bioinformatics/btz305 [↩](#cite-36)
37. <a id="ref-37"></a> Togkousidis, Anastasis, Oleksiy M. Kozlov, Julia Haag, Dimitri Höhler, and Alexandros Stamatakis. 2023. "Adaptive RAxML-NG: Accelerating phylogenetic inference under maximum likelihood using dataset difficulty." _Molecular Biology and Evolution_ 40 (10): msad227. https://doi.org/10.1093/molbev/msad227 [↩](#cite-37)
38. <a id="ref-38"></a> Haag, Julia, Dimitri Höhler, Ben Bettisworth, and Alexandros Stamatakis. 2022. "From easy to hopeless: Predicting the difficulty of phylogenetic analyses." _Molecular Biology and Evolution_ 39 (12): msac254. https://doi.org/10.1093/molbev/msac254 [↩¹](#cite-38a) [↩²](#cite-38b)
39. <a id="ref-39"></a> Kozlov, Oleksiy M., Anastasis Togkousidis, Christoph Stelz, Dimitri Höhler, Julius Wiegert, and Alexandros Stamatakis. 2026. "RAxML-NG 2: Automatic model selection, novel tree search heuristics, and fast branch support metrics." Preprint. https://doi.org/10.64898/2026.09.09.750097 [↩](#cite-39)
40. <a id="ref-40"></a> Price, Morgan N., Paramvir S. Dehal, and Adam P. Arkin. 2010. "FastTree 2: Approximately maximum-likelihood trees for large alignments." _PLoS ONE_ 5 (3): e9490. https://doi.org/10.1371/journal.pone.0009490 [↩](#cite-40)
41. <a id="ref-41"></a> Piñeiro, César, and Juan C. Pichel. 2024. "Efficient phylogenetic tree inference for massive taxonomic datasets: Harnessing the power of a server to analyze 1 million taxa." _GigaScience_ 13: giae055. https://doi.org/10.1093/gigascience/giae055 [↩](#cite-41)
42. <a id="ref-42"></a> Guindon, Stéphane, Jean-François Dufayard, Vincent Lefort, Maria Anisimova, Wim Hordijk, and Olivier Gascuel. 2010. "New algorithms and methods to estimate maximum-likelihood phylogenies: Assessing the performance of PhyML 3.0." _Systematic Biology_ 59 (3): 307-321. https://doi.org/10.1093/sysbio/syq010 [↩](#cite-42)
43. <a id="ref-43"></a> Zhou, Xiaofan, Xing-Xing Shen, Chris Todd Hittinger, and Antonis Rokas. 2018. "Evaluating fast maximum likelihood-based phylogenetic programs using empirical phylogenomic data sets." _Molecular Biology and Evolution_ 35 (2): 486-503. https://doi.org/10.1093/molbev/msx302 [↩](#cite-43)
44. <a id="ref-44"></a> Whelan, Simon, and Daniel Money. 2010. "The prevalence of multifurcations in tree-space and their implications for tree-search." _Molecular Biology and Evolution_ 27 (12): 2674-2677. https://doi.org/10.1093/molbev/msq163 [↩](#cite-44)
45. <a id="ref-45"></a> Kalyaanamoorthy, Subha, Bui Quang Minh, Thomas K. F. Wong, Arndt von Haeseler, and Lars S. Jermiin. 2017. "ModelFinder: Fast model selection for accurate phylogenetic estimates." _Nature Methods_ 14 (6): 587-589. https://doi.org/10.1038/nmeth.4285 [↩](#cite-45)
46. <a id="ref-46"></a> Abadi, Shiran, Dana Azouri, Tal Pupko, and Itay Mayrose. 2019. "Model selection may not be a mandatory step for phylogeny reconstruction." _Nature Communications_ 10: 934. https://doi.org/10.1038/s41467-019-08822-w [↩](#cite-46)
47. <a id="ref-47"></a> De Maio, Nicola, Myrthe Willemsen, Samuel Martin, et al. 2026. "Rate variation and recurrent sequence errors in pandemic-scale phylogenetics." _Nature Methods_ 23 (3): 565-573. https://doi.org/10.1038/s41592-025-02932-8 [↩](#cite-47)
48. <a id="ref-48"></a> Ly-Trong, Nhan, Samuel Martin, Nick Goldman, Nicola De Maio, and Bui Quang Minh. 2026. "CMAPLE 2: Fast and accurate phylogenetic inference for millions of pathogen genomes." Preprint. https://doi.org/10.64898/2026.06.15.732229 [↩](#cite-48)
49. <a id="ref-49"></a> Matsen, Frederick A., Robin B. Kodner, and E. Virginia Armbrust. 2010. "pplacer: Linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree." _BMC Bioinformatics_ 11: 538. https://doi.org/10.1186/1471-2105-11-538 [↩](#cite-49)
50. <a id="ref-50"></a> Barbera, Pierre, Alexey M. Kozlov, Lucas Czech, et al. 2019. "EPA-ng: Massively parallel evolutionary placement of genetic sequences." _Systematic Biology_ 68 (2): 365-369. https://doi.org/10.1093/sysbio/syy054 [↩](#cite-50)
51. <a id="ref-51"></a> Balaban, Metin, Shahab Sarmashghi, and Siavash Mirarab. 2020. "APPLES: Scalable distance-based phylogenetic placement with or without alignments." _Systematic Biology_ 69 (3): 566-578. https://doi.org/10.1093/sysbio/syz063 [↩](#cite-51)
52. <a id="ref-52"></a> Wedell, Eleanor, Yirong Cai, and Tandy Warnow. 2023. "SCAMPP: Scaling alignment-based phylogenetic placement to large trees." _IEEE/ACM Transactions on Computational Biology and Bioinformatics_ 20 (2): 1417-1430. https://doi.org/10.1109/tcbb.2022.3170386 [↩](#cite-52)
53. <a id="ref-53"></a> Felsenstein, Joseph. 1985. "Confidence limits on phylogenies: An approach using the bootstrap." _Evolution_ 39 (4): 783-791. https://doi.org/10.1111/j.1558-5646.1985.tb00420.x [↩](#cite-53)
54. <a id="ref-54"></a> Hoang, Diep Thi, Olga Chernomor, Arndt von Haeseler, Bui Quang Minh, and Le Sy Vinh. 2018. "UFBoot2: Improving the ultrafast bootstrap approximation." _Molecular Biology and Evolution_ 35 (2): 518-522. https://doi.org/10.1093/molbev/msx281 [↩](#cite-54)
55. <a id="ref-55"></a> Anisimova, Maria, and Olivier Gascuel. 2006. "Approximate likelihood-ratio test for branches: A fast, accurate, and powerful alternative." _Systematic Biology_ 55 (4): 539-552. https://doi.org/10.1080/10635150600755453 [↩](#cite-55)
56. <a id="ref-56"></a> Anisimova, Maria, Manuel Gil, Jean-François Dufayard, Christophe Dessimoz, and Olivier Gascuel. 2011. "Survey of branch support methods demonstrates accuracy, power, and robustness of fast likelihood-based approximation schemes." _Systematic Biology_ 60 (5): 685-699. https://doi.org/10.1093/sysbio/syr041 [↩](#cite-56)
57. <a id="ref-57"></a> Lemoine, F., J.-B. Domelevo Entfellner, E. Wilkinson, et al. 2018. "Renewing Felsenstein's phylogenetic bootstrap in the era of big data." _Nature_ 556 (7702): 452-456. https://doi.org/10.1038/s41586-018-0043-0 [↩](#cite-57)
58. <a id="ref-58"></a> Baele, Guy, Xiang Ji, Gabriel W. Hassler, et al. 2025. "BEAST X for Bayesian phylogenetic, phylogeographic and phylodynamic inference." _Nature Methods_ 22 (8): 1653-1656. https://doi.org/10.1038/s41592-025-02751-x [↩](#cite-58)
59. <a id="ref-59"></a> Bouckaert, Remco R., Paula H. Weidemüller, Luis R. Esquivel Gomez, and Nicola F. Müller. 2025. "Improving the scalability of Bayesian phylodynamic inference through efficient MCMC proposals." _bioRxiv_ preprint. https://doi.org/10.1101/2025.06.18.660471 [↩](#cite-59)
60. <a id="ref-60"></a> Gill, Mandev S., Philippe Lemey, Marc A. Suchard, Andrew Rambaut, and Guy Baele. 2020. "Online Bayesian phylodynamic inference in BEAST with application to epidemic reconstruction." _Molecular Biology and Evolution_ 37 (6): 1832-1842. https://doi.org/10.1093/molbev/msaa047 [↩](#cite-60)
61. <a id="ref-61"></a> Gao, Jiansi, Marius Brusselmans, Luiz M. Carvalho, Marc A. Suchard, Guy Baele, and Frederick A. Matsen. 2026. "Biological causes and impacts of rugged tree landscapes in phylodynamic inference." _Proceedings of the National Academy of Sciences_ 123 (2): e2510938123. https://doi.org/10.1073/pnas.2510938123 [↩](#cite-61)
62. <a id="ref-62"></a> To, Thu-Hien, Matthieu Jung, Samantha Lycett, and Olivier Gascuel. 2016. "Fast dating using least-squares criteria and algorithms." _Systematic Biology_ 65 (1): 82-97. https://doi.org/10.1093/sysbio/syv068 [↩](#cite-62)
63. <a id="ref-63"></a> Sanderson, Theo. 2021. "Chronumental: Time tree estimation from very large phylogenies." _bioRxiv_ preprint. https://doi.org/10.1101/2021.10.27.465994 [↩](#cite-63)
64. <a id="ref-64"></a> Ji, Xiang, Zhenyu Zhang, Andrew Holbrook, et al. 2020. "Gradients do grow on trees: A linear-time O(N)-dimensional gradient for statistical phylogenetics." _Molecular Biology and Evolution_ 37 (10): 3047-3060. https://doi.org/10.1093/molbev/msaa130 [↩](#cite-64)
65. <a id="ref-65"></a> Turakhia, Yatish, Bryan Thornlow, Angie Hinrichs, et al. 2022. "Pandemic-scale phylogenomics reveals the SARS-CoV-2 recombination landscape." _Nature_ 609 (7929): 994-997. https://doi.org/10.1038/s41586-022-05189-9 [↩](#cite-65)
66. <a id="ref-66"></a> Zhan, Shing H., Yan Wong, Anastasia Ignatieva, et al. 2023. "A pandemic-scale ancestral recombination graph for SARS-CoV-2." _bioRxiv_ preprint. https://doi.org/10.1101/2023.06.08.544212 [↩](#cite-66)
67. <a id="ref-67"></a> Kelleher, Jerome, Yan Wong, Anthony W. Wohns, Chaimaa Fadil, Patrick K. Albers, and Gil McVean. 2019. "Inferring whole-genome histories in large population datasets." _Nature Genetics_ 51 (9): 1330-1338. https://doi.org/10.1038/s41588-019-0483-y [↩](#cite-67)
68. <a id="ref-68"></a> Croucher, Nicholas J., Andrew J. Page, Thomas R. Connor, et al. 2015. "Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins." _Nucleic Acids Research_ 43 (3): e15. https://doi.org/10.1093/nar/gku1196 [↩](#cite-68)
69. <a id="ref-69"></a> Didelot, Xavier, and Daniel J. Wilson. 2015. "ClonalFrameML: Efficient inference of recombination in whole bacterial genomes." _PLOS Computational Biology_ 11 (2): e1004041. https://doi.org/10.1371/journal.pcbi.1004041 [↩](#cite-69)
70. <a id="ref-70"></a> Nesterenko, Luca, Luc Blassel, Philippe Veber, Bastien Boussau, and Laurent Jacob. 2025. "Phyloformer: Fast, accurate, and versatile phylogenetic reconstruction with deep neural networks." _Molecular Biology and Evolution_ 42 (4): msaf051. https://doi.org/10.1093/molbev/msaf051 [↩](#cite-70)
71. <a id="ref-71"></a> Zhang, Xinru, Shizhe Ding, Chungong Yu, Jianquan Zhao, and Dongbo Bu. 2025. "Accurate and efficient phylogenetic inference through end-to-end deep learning." _Molecular Biology and Evolution_ 42 (11): msaf260. https://doi.org/10.1093/molbev/msaf260 [↩](#cite-71)
72. <a id="ref-72"></a> Yang, Chen, Zixin Zhuang, Piyumal Demotte, et al. 2027. "IQ-NET: Fast and accurate quartet phylogenetic inference using deep learning trained on empirical DNA alignments." _Molecular Phylogenetics and Evolution_ 226: 108744. https://doi.org/10.1016/j.ympev.2026.108744 [↩](#cite-72)
73. <a id="ref-73"></a> Azouri, Dana, Shiran Abadi, Yishay Mansour, Itay Mayrose, and Tal Pupko. 2021. "Harnessing machine learning to guide phylogenetic-tree search algorithms." _Nature Communications_ 12: 1983. https://doi.org/10.1038/s41467-021-22073-8 [↩](#cite-73)
74. <a id="ref-74"></a> Azouri, Dana, Oz Granit, Michael Alburquerque, Yishay Mansour, Tal Pupko, and Itay Mayrose. 2024. "The tree reconstruction game: Phylogenetic reconstruction using reinforcement learning." _Molecular Biology and Evolution_ 41 (6): msae105. https://doi.org/10.1093/molbev/msae105 [↩](#cite-74)
75. <a id="ref-75"></a> Fourment, Mathieu, Christiaan J. Swanepoel, Jared G. Galloway, et al. 2023. "Automatic differentiation is no panacea for phylogenetic gradient computation." _Genome Biology and Evolution_ 15 (6): evad099. https://doi.org/10.1093/gbe/evad099 [↩](#cite-75)
76. <a id="ref-76"></a> Thornlow, Bryan, Alexander Kramer, Cheng Ye, et al. 2021. "Online phylogenetics using parsimony produces slightly better trees and is dramatically more efficient for large SARS-CoV-2 phylogenies than de novo and maximum-likelihood approaches." _bioRxiv_ preprint. https://doi.org/10.1101/2021.12.02.471004 [↩](#cite-76)
