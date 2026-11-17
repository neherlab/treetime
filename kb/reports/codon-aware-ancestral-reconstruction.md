# Codon-aware ancestral reconstruction

This report collects the evidence for one goal: TreeTime v1 must write ancestral nucleotide and amino-acid sequences that agree, so that the translation of every node's nucleotide sequence equals that node's amino-acid sequence. It covers the two ways in which Nextstrain pipelines get ancestral amino acids, the disagreement of the two ways measured in v1, the theory of turning nucleotide posteriors into codons, what established reconstruction tools do, annotations and overlapping reading frames, the evidence on codon substitution models, the cost of the options, and pitfalls found in an unmerged prototype.

The design options and the recommended design are in [kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md).

Code facts describe the `rust` branch as of 2026-11-12. Links to third-party code are pinned to the revision that was read. TreeTime v0 and augur are legacy implementations here: the report describes their behavior where Nextstrain consumers depend on it, and parity with them is no argument for a design.

## Summary

- **Consistent outputs come from decoding, not from a codon model.** When each node gets one codon per CDS codon, and the nucleotide and amino-acid outputs both come from that codon, the outputs agree by construction. A codon model alone does not give this: a decoder that writes the most probable codon as nucleotides and the most probable amino acid as the residue disagrees with itself (section [Theory of codon decoding](#theory-of-codon-decoding))
- **The codon posterior is already available.** Under site-independent nucleotide models, which v1 uses today, the joint posterior of the three nucleotides of a codon at a node is exactly the product of the three site posteriors
- **The two routes almost agree in v1 today**: 76 of 9.1 million internal-node residues differ with both sides determined on `rsv/a/2000`, and 17 of 110 million on `mpox/clade-ii/2000`. Model differences at confident codons and gap handling cause these differences, not the order of reconstruction and translation
- **Real ancestral stop codons exist**: SARS-CoV-2 ORF8 Q27* sits on an internal node of clade 20I, and 57 internal nodes of `rsv/a/2000` and 142 of `mpox/clade-ii/2000` carry stop codons that their tips share. A rule that forbids stop codons, and a 61-state codon model, overwrite them
- **Overlapping reading frames are the normal case**: 46 of 120 virus annotations have them, among them every annotated dataset in `data/`. No published method reconstructs ancestral states across overlapping frames
- **The annotation reader does not read the Nextstrain annotations correctly**: it rejects strand `.` (all RSV rows), reads no file with `gene` rows only (SARS-CoV-2), ignores phase and `Parent`, and does not check that alignment columns are reference coordinates
- **No study shows that codon models give better ancestral states at viral divergence**: the one direct benchmark is a tie, 98.4% (nucleotide model) against 98.3% (codon model) per-site accuracy <a id="cite-1a"></a>[Randall et al. 2016](https://doi.org/10.1038/ncomms12847) [[1](#ref-1)]
- **The amino-acid route is expensive**: with one job it adds 3.0 times the run time of the nucleotide reconstruction on `mpox/clade-ii/2000`, and 2.7 times on `rsv/a/2000`
- **`--aa-model jtt92` writes invented residues** at the annotated stop codon of every internal node, where every tip has `X`

## Background

### Two routes to ancestral amino acids in Nextstrain

Pipelines get ancestral amino acids in one of two ways. This report calls them the translation route and the peptide route.

- **Translation route**: `augur translate` translates the nucleotide reconstruction of every node. 10 of the 13 Nextstrain pipeline repositories that were checked use it for their main build
- **Peptide route**: `augur ancestral --translations` runs a second, independent reconstruction in amino-acid space, from tip peptides that Nextclade translated. seasonal-flu and ncov use it (ncov with its own script), and ebola uses it for two of its three species
- **Both routes**: the RSV pipeline runs both, and the output of the translation route overwrites the output of the peptide route at export ([export.smk#L4-L5](https://github.com/nextstrain/rsv/blob/7d91d8fd0ff98d3a33ad8191a8ad712d05904e58/workflow/snakemake_rules/export.smk#L4-L5))

`augur ancestral` seeds each reconstruction with JC69 and then infers a GTR model (`infer_gtr=True`) for both alphabets ([ancestral.py#L302-L307](https://github.com/nextstrain/augur/blob/0dc610704848c3a472f7715a1161277ef9b625a4/augur/ancestral.py#L302-L307)). Only the ncov script uses JC69 without inference.

### Why the peptide route disagrees with the nucleotide reconstruction

The two reconstructions use different models, alphabets and decoding rules. In addition, Nextclade peptides do not line up codon by codon with the nucleotide alignment: Nextclade strips gaps, aligns each peptide to the reference peptide, strips insertions, and masks frameshifted codons with `X`. Two consequences follow:

- A node can carry an amino acid that its own codon does not code
- A branch can carry an amino-acid mutation without any nucleotide change in that codon

### What Auspice expects

Auspice derives the ancestral state of a position from the "from" letter of the first mutation that it meets on the path from the root, and checks it against the root sequence. A mismatch only logs `console.error` ([setGenotype.js](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/util/setGenotype.js#L3)). Mutation chains that do not match the parent states therefore give wrong colors without a visible error.

### TreeTime v1 today

- **TreeTime implements the peptide route only.** `treetime ancestral --translations <peptides> --annotation <gff3>` runs one independent amino-acid partition per CDS (`fn reconstruct_aa()` [packages/treetime/src/ancestral/aa.rs#L25](../../packages/treetime/src/ancestral/aa.rs#L25)), each with its own GTR model and fitted rate ([packages/treetime/src/ancestral/aa.rs#L39-L68](../../packages/treetime/src/ancestral/aa.rs#L39-L68)). This path ignores `--method-anc` and `--include-leaves` ([kb/issues/M-ancestral-aa-path-ignores-method-and-include-leaves.md](../issues/M-ancestral-aa-path-ignores-method-and-include-leaves.md))
- **No codon substitution model exists**
- **Partitions are site sets of the `--alignment` files**, including codon-position partitions such as `a-b\3` ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)). `--translations` and `--annotation` accept one alignment file
- **Annotated datasets**: `data/rsv/a/*`, `data/mpox/clade-ii/*` and `data/sc2/*` have an `annotation.gff3`. Peptide files exist for `rsv/a/20`, `rsv/a/2000` and `mpox/clade-ii/2000`, made by the `data/translate` script

## Measured disagreement of the two routes

### Method

`treetime ancestral` ran on the `rust` branch of 2026-11-12, release build, with default settings: marginal reconstruction, inferred GTR, sparse representation, `--gap-fill only-terminal`. Each dataset ran once without and once with `--translations`. The nucleotide output and the tree are byte-identical in both runs, so both routes work on the same tree. The CDS of every internal node was translated from the nucleotide output and compared residue by residue with the amino-acid output.

- **Translation rule of the comparison**: `---` gives `-`; any other codon with a gap gives `X`; an ambiguous codon gives its amino acid when all resolutions code the same amino acid, else `X`. The annotated terminal stop codon is counted apart
- **Peptide input**: the bundled peptide files are translations of the same alignment by `data/translate`, so peptide columns match the CDS codons exactly, and the comparison contains no column misalignment from indels
- **Annotations**: the alignment length equals the reference length (15,225 and 197,209 nucleotides), every CDS length is a multiple of 3, and each CDS is one GFF3 row with phase 0. The RSV runs used `--cds-selection`, because the RSV annotation used strand `.`. The mpox runs used `--annotation`: 179 CDS, 96 of them on the minus strand

Commands, with `D=data/<dataset>` and `O=<output folder>`:

```
treetime ancestral --tree $D/tree.nwk --alignment $D/aln.fasta.xz --output-sequences $O/nuc.nuc.fasta --output-tree-nwk $O/nuc.tree.nwk
treetime ancestral --tree $D/tree.nwk --alignment $D/aln.fasta.xz <CDS> --translations "$D/translations/{cds}.fasta.xz" --output-sequences $O/infer.nuc.fasta --output-translations "$O/infer.aa/{cds}.fasta" --output-tree-nwk $O/infer.tree.nwk
treetime ancestral --tree $D/tree.nwk --alignment $D/aln.fasta.xz --ambiguity-threshold 0.5 --output-sequences $O/amb05.nuc.fasta

<CDS> = --cds-selection NS1,NS2,N,P,M,SH,G,F,M2-1,M2-2,L   (RSV)
<CDS> = --annotation $D/annotation.gff3                   (mpox)
```

### Results

| Dataset | Internal-node residues | Both determined, different | One side `X` or `-`, different | Internal stop codons (node-codons, both routes) |
| --- | --- | --- | --- | --- |
| `rsv/a/20` | 82,026 | 11 (0.013%) | 3 | 0 |
| `rsv/a/2000` | 9,100,329 | 76 (0.0008%) | 496 (0.0055%) | 57 (3 columns, all in G) |
| `mpox/clade-ii/2000` | 109,640,673 | 17 (0.00002%) | 4,915 (0.0045%) | 142 (13 columns) |

- **Model differences, not decoding, cause most both-determined differences.** 55 of the 76 RSV differences and 15 of the 17 mpox differences sit at codons that the run with `--ambiguity-threshold 0.5` does not flag: at each of their positions the most probable nucleotide leads every other nucleotide by more than 0.5. One RSV column (P codon 176, D against E, a difference at the third position) holds 34 of the 76. The nucleotide route fits one GTR model to the genome, the amino-acid route one per CDS, so the two are different estimators
- **Gap handling causes most other differences.** `--gap-fill only-terminal` turns the leading and trailing gaps of each sequence into unknown characters. The amino-acid route applies it to each peptide (`fn read_cds_translations()` [packages/app-commands/src/commands/ancestral/run.rs#L376-L397](../../packages/app-commands/src/commands/ancestral/run.rs#L376-L397)), so a deletion that touches the start or end of a CDS becomes `X` there, while it stays an interior gap of the genome in the nucleotide route. This explains most of the 4,915 mpox differences: the CDS with many differences of this kind all have tips with such gap runs (11,413 tip codons on mpox)
- **Internal stop codons are inherited.** Every internal stop codon sits in a column where tips of the same clade carry the same stop codon, and both routes reconstruct the same 57 and 142 internal stop codons
- **No order effect is visible.** Only 16 internal nodes of `rsv/a/2000` and 1 of `mpox/clade-ii/2000` carry a translated residue that none of their determined children carries
- **`--aa-model jtt92` invents residues.** Its 20-state alphabet has no stop symbol (`AlphabetName::AaNoStop`, [packages/app-commands/src/commands/ancestral/aa_model.rs#L30-L34](../../packages/app-commands/src/commands/ancestral/aa_model.rs#L30-L34)). Every tip has `X` at the annotated stop codon, but every internal node gets `M` there in F, OPG001 and OPG105, and `Q` in L: 1,997 to 1,999 nodes per CDS on RSV and 1,987 on mpox. With `jtt92` the route also reconstructs none of the internal stop codons

> [!IMPORTANT]
> **Investigation required.** RSV L has 136 differences of the gap class but only 10 tip codons with CDS-terminal gaps, so the terminal-gap fill does not explain them. Their cause is not traced. The mechanism by which `jtt92` writes `M` or `Q` at the stop codon is also not traced.

> [!NOTE]
> Both trees are shallow outbreak trees. Deeper datasets (influenza, dengue, tb) can have more diffuse codons and more model-driven differences. They have no annotation in `data/`, so they were not measured.

### Share of uncertain codons

No CLI output contains the site posteriors, so the share of <a id="gloss-use-1"></a>diffuse <sup>[1](#gloss-1)</sup> codons, where the decoding rules can differ, was not measured directly. `--ambiguity-threshold t` reports every state within `t` of the most probable state ([packages/treetime/src/partition/marginal/posterior_states.rs#L34-L52](../../packages/treetime/src/partition/marginal/posterior_states.rs#L34-L52)). With `t = 0.5` it flags every position whose most probable nucleotide leads the second by at most 0.5, which includes every position with a posterior below 1/2.

| Dataset | Internal-node codon sites | Codon sites with a flagged position | Columns with a flagged position |
| --- | --- | --- | --- |
| `rsv/a/2000` | 9,122,296 | 448 (0.0049%) | 250 of 4,568 |
| `mpox/clade-ii/2000` | 109,996,346 | 27 (0.000025%) | 16 of 55,358 |

These counts include every codon site that has a position below 1/2. They do not bound the share of diffuse codons, because a codon can be diffuse when each position is above the threshold: three positions at 0.7 give a most probable triplet of 0.34. They show the order of magnitude on shallow trees.

> [!IMPORTANT]
> **Investigation required.** Measure the share of diffuse codons, and the share where the decoding rules choose different codons, with a decoder that reads the site posteriors, on the annotated datasets and on a deep tree.

### Cost of the amino-acid route

Means of 3 runs, release build, spread under 3%:

| Dataset | Run | Wall time, 20 jobs | CPU time (user), 20 jobs | Wall time, 1 job | Peak memory |
| --- | --- | --- | --- | --- | --- |
| `rsv/a/2000` | nucleotides only | 0.32 s | 1.79 s | 1.08 s | 276 MB |
| `rsv/a/2000` | with amino acids, inferred GTR | 1.02 s | 7.89 s | 3.95 s | 314 MB |
| `rsv/a/2000` | with amino acids, `jtt92` | 0.98 s | 7.52 s | 3.73 s | 307 MB |
| `mpox/clade-ii/2000` | nucleotides only | 3.87 s | 29.1 s | 14.5 s | 2.59 GB |
| `mpox/clade-ii/2000` | with amino acids, inferred GTR | 15.0 s | 126.5 s | 57.9 s | 2.61 GB |
| `mpox/clade-ii/2000` | with amino acids, `jtt92` | 14.6 s | 121.2 s | 55.3 s | 2.59 GB |

- The amino-acid route adds 2.7 times (RSV) and 3.0 times (mpox) the time of the nucleotide run with one job
- `jtt92` costs only 4% less than the inferred GTR, so model fitting is not the main cost. The CDS loop runs one CDS after the other ([packages/treetime/src/ancestral/aa.rs#L39](../../packages/treetime/src/ancestral/aa.rs#L39))
- Peak memory comes from the nucleotide part; the amino-acid route adds at most about 40 MB

## Theory of codon decoding

### Notation

For one CDS codon at one node:

- $p(c)$: posterior of triplet $c$, over the 64 triplets, or over the 61 <a id="gloss-use-2"></a>sense <sup>[2](#gloss-2)</sup> codons under a codon model
- $g$: the genetic code, which maps a triplet to an amino acid or to the stop symbol `*`
- $q(a) = \sum_{g(c)=a} p(c)$: posterior of amino acid $a$, or of the stop symbol
- $r_k(n)$: posterior of nucleotide $n$ at codon position $k \in \lbrace 1, 2, 3 \rbrace$
- $\hat c$: the most probable triplet
- <a id="gloss-use-3"></a>per-site argmax <sup>[3](#gloss-3)</sup>: the triplet of the most probable nucleotide at each position, which TreeTime writes today

Decoding is the step that turns the codon posterior of one node into one codon. When the nucleotide output and the amino-acid output both come from that one codon, they agree at every node, and every amino-acid change on a branch comes with a nucleotide change in that codon on the same branch, because translation is a function.

### The codon posterior of a site-independent model

When sites evolve independently given the tree, the branch lengths and the model parameters, the posterior of all ancestral states factorizes over sites. The joint posterior of the three nucleotides of a codon at one node is then the product of the three site posteriors:

$$p(c) = r_1(c_1)\, r_2(c_2)\, r_3(c_3)$$

- **Checked numerically**: brute-force enumeration of all 64 triplets on a 4-tip tree matches the product to within 5.4e-15
- **Holds** for one GTR model, for codon-position partitions (`a-b\3`), and for rate categories that are summed per site
- **Fails** for rates that are correlated between neighboring sites, for context-dependent and overlapping-frame models <a id="cite-2a"></a>[Pedersen and Jensen 2001](https://doi.org/10.1093/oxfordjournals.molbev.a003859) [[2](#ref-2)], under full Bayesian integration over tree and parameters, and for codon models
- **Holds in v1 today**: v1 has fixed per-site rates and no rate categories (`site_rates` in [packages/treetime/src/gtr/gtr.rs](../../packages/treetime/src/gtr/gtr.rs))

PAML builds amino-acid posteriors in the same way for coding data that its nucleotide program analyzes: it multiplies the three site posteriors of a codon, drops the stop mass and renormalizes. Its source comment calls the result "approximate" because the stop codons are ignored ([treesub.c#L6389-L6415](https://github.com/abacus-gene/paml/blob/4c7902fe972737ef5e80bb18f159a6e6acace3d3/src/treesub.c#L6389-L6415)).

### Decoding rules from a loss function

Each decoding rule minimizes an expected loss. Take the loss $L = \lambda \cdot \mathbb{1}[\text{amino acid wrong}] + \text{Hamming}(\text{nucleotides})$. Its <a id="gloss-use-4"></a>Bayes rule <sup>[4](#gloss-4)</sup> picks the triplet that maximizes

$$\lambda\, q(g(c)) + \sum_{k=1}^{3} r_k(c_k)$$

with the stop symbol as a 21st value of $g$. Four rules follow:

| Rule | Loss | Choice |
| --- | --- | --- |
| Amino acid first | $\lambda \to \infty$ (<a id="gloss-use-5"></a>lexicographic loss <sup>[5](#gloss-5)</sup>) | the most probable amino acid or stop, then its codon with the largest $\sum_k r_k$ |
| Most probable codon | codon 0-1 loss | $\hat c$ |
| Nucleotide first | $\lambda = 0$, sense codons only (the constrained <a id="gloss-use-6"></a>centroid estimator <sup>[6](#gloss-6)</sup> <a id="cite-3a"></a>[Carvalho and Lawrence 2008](https://doi.org/10.1073/pnas.0712329105) [[3](#ref-3)]) | the sense codon with the largest $\sum_k r_k$ |
| Per-site argmax | $\lambda = 0$, all 64 triplets | the most probable nucleotide at each position, translated |

The amino-acid-first rule is already exact for every finite $\lambda > 3 / \Delta q$, where $\Delta q$ is the gap between the two largest values of $q$, because $\sum_k r_k$ lies between 0 and 3. No paper states this derivation for codon reconstruction.

A fifth rule, the split decoder, takes the nucleotides from $\hat c$ and the amino acid from the largest $q$. It is the only rule here that gives inconsistent outputs: with CTG 0.40 (Leu), ATT 0.35 and ATC 0.25 (Ile), the nucleotides say Leu and the amino acid says Ile (0.60).

### When the rules agree

If $p(\hat c) \ge 1/2$, the amino acid of $\hat c$ has $q \ge 1/2$, and each nucleotide of $\hat c$ has a site posterior of at least 1/2. All four rules then choose $\hat c$, apart from ties at exactly 1/2. The rules can differ only at diffuse codons, where $p(\hat c) < 1/2$ and a competing amino acid spreads its mass over at least two synonymous codons.

Under the product-form posterior of a site-independent model, the per-site argmax is also the most probable triplet, and it has the largest $\sum_k r_k$ of all 64 triplets. So the most-probable-codon rule and the nucleotide-first rule reproduce today's nucleotide output, except where that triplet is a stop codon. The amino-acid-first rule keeps today's output whenever the per-site argmax codes the most probable amino acid, and changes nucleotides only where it codes a less probable one. A check over 60,000 random cases found no counterexample.

### Worked example under a product-form posterior

Position 1 is A (0.55) or C (0.45), position 2 is T (1.0), position 3 is G (0.50), A (0.20), C (0.15) or T (0.15).

- **Per-site argmax**: ATG (Met), with $p = q(\text{Met}) = 0.275$. Today's output, translated as the translation route does, reports Met
- **Most probable amino acid**: Leu, with $q = 0.45$ spread over the four CTN codons
- **Amino acid first**: CTG (Leu). Expected amino-acid error 0.55 against 0.725 for ATG; expected nucleotide error 1.05 against 0.95

The translated per-site argmax reports an amino acid that is clearly less probable than another one. The amino-acid-first rule removes that error at the cost of 0.1 expected nucleotide errors.

### Worked example under a codon model

Posteriors under a codon model need not be products of site posteriors. Example: TAC (Tyr) 0.30, TTA (Leu) 0.25, CTA (Leu) 0.20, CAA (Gln) 0.25.

| Rule | Codon | Expected amino-acid error | Expected nucleotide error |
| --- | --- | --- | --- |
| Most probable codon | TAC (Tyr) | 0.70 | 1.60 |
| Amino acid first | TTA (Leu) | 0.55 | 1.30 |
| Nucleotide first | TTA (Leu) or CAA (Gln), tied | 0.55 or 0.75 | 1.30 |
| Per-site argmax | TAA (stop) | 1.00 | 1.20 |

The site posteriors are T 0.55 at position 1, A 0.55 at position 2 and A 0.70 at position 3, so the per-site argmax gives the stop codon TAA, which has no posterior mass at all.

### Stop codons

- **Renormalizing over sense codons does not change the amino-acid decision.** It divides the $q$ of every sense amino acid by the same $1 - q(\ast)$, so their order stays. It changes only the reported probabilities and, when renormalized site posteriors feed the codon choice, the codon (23 to 497 of 20,000 random cases)
- **Real ancestral stop codons exist.** In the Nextclade SARS-CoV-2 reference tree, the internal node `NODE_0001012` of clade 20I, with 5 children, carries `ORF8: Q27*` ([tree.json](https://github.com/nextstrain/nextclade_data/blob/a231e7badf4747b9665708f7f662079b6f94264e/data/nextstrain/sars-cov-2/wuhan-hu-1/orfs/tree.json)). A rule that forbids stop codons at internal nodes forces a sense codon there, and spreads Q27* as repeated changes over the five child branches
- **Exact conditioning on "no internal stop codon"** is possible with 64-state pruning: the transition matrix is the Kronecker product $P_1(t) \otimes P_2(t) \otimes P_3(t)$ of the nucleotide model, with a stop indicator at internal nodes. It changes nodes next to stop-prone nodes as well (codon changes at 5 to 34 of 300 nodes on constructed trees), costs about 16 times more per affected codon (768 against 48 multiply-adds per edge), and forbids the real stop codons above

### Histories

Per-node decoding gives a consistent pair of outputs at each node, but it does not give a parsimonious history: it can produce a residue that no child carries, or extra synonymous changes, where a joint history would differ. Today's per-site marginal reconstruction has the same property. Joint and marginal reconstructions do differ in practice: 13 differences at 5 positions over 19 internal nodes of 21 cytochrome b sequences <a id="cite-4a"></a>[Pupko et al. 2000](https://doi.org/10.1093/oxfordjournals.molbev.a026369) [[4](#ref-4)]. The <a id="gloss-use-7"></a>joint reconstruction <sup>[7](#gloss-7)</sup> of codons needs a joint pass that v1 does not have. The most probable amino-acid history is a <a id="gloss-use-8"></a>marginal MAP <sup>[8](#gloss-8)</sup> problem; whether the tree structure allows an exact algorithm is not known.

Point estimates also bias counts: single best reconstructions bias composition and substitution counts <a id="cite-5"></a>[Williams et al. 2006](https://doi.org/10.1371/journal.pcbi.0020069) [[5](#ref-5)]; <a id="cite-6"></a>[Krishnan et al. 2004](https://doi.org/10.1093/molbev/msh198) [[6](#ref-6)]; <a id="cite-7"></a>[Matsumoto et al. 2015](https://doi.org/10.1534/genetics.115.177386) [[7](#ref-7)]. Expected counts from edge posteriors or <a id="gloss-use-9"></a>robust counting <sup>[9](#gloss-9)</sup> <a id="cite-8a"></a>[Minin and Suchard 2008](https://doi.org/10.1007/s00285-007-0120-8) [[8](#ref-8)]; <a id="cite-9a"></a>[Lemey et al. 2012](https://doi.org/10.1093/bioinformatics/bts580) [[9](#ref-9)] are the sound basis for synonymous and nonsynonymous summaries. Codon-level minimum posterior expected error sets exist for uncertainty reporting <a id="cite-10"></a>[Oliva et al. 2019](https://doi.org/10.1093/bioinformatics/btz249) [[10](#ref-10)], but they need a loss over codon sets that no source defines.

## How established tools reconstruct and translate

| Tool | Ancestral coding states | Source |
| --- | --- | --- |
| PAML | Prints the most probable codon and, beside it, the separately summed most probable amino acid, and chooses neither | [treesub.c#L6537](https://github.com/abacus-gene/paml/blob/4c7902fe972737ef5e80bb18f159a6e6acace3d3/src/treesub.c#L6537) |
| HyPhy | Translates joint maximum-likelihood codons | [selection_lib.ibf#L68](https://github.com/veg/hyphy/blob/1fa8fd0469328b8d6df2c7273077170765beee5f/res/TemplateBatchFiles/SelectionAnalyses/modules/selection_lib.ibf#L68) |
| IQ-TREE 3 | Marginal codon states only, no amino-acid output; a state becomes unknown below `--asr-min` | [treetesting.cpp#L230-L330](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/main/treetesting.cpp#L230-L330), [phylotreesse.cpp#L1488-L1509](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/tree/phylotreesse.cpp#L1488-L1509) |
| BEAST 1 | Draws codons from the joint posterior by <a id="gloss-use-10"></a>pre-order sampling <sup>[10](#gloss-10)</sup> | [AncestralStateBeagleTreeLikelihood.java#L414-L530](https://github.com/beast-dev/beast-mcmc/blob/d1a0c7ef52ea1987c869e0b34889ffac62f4da8a/src/dr/evomodel/treelikelihood/AncestralStateBeagleTreeLikelihood.java#L414-L530) |
| UShER `matUtils translate` | Derives amino acids from the nucleotide states of the tree | [translate.cpp#L240-L294](https://github.com/yatisht/usher/blob/ac9c982d937c3bc43e2c16a0b73d30cf2b937118/src/matUtils/translate.cpp#L240-L294) |

No tool enforces consistent outputs by a rule other than translating one codon.

### Stop codons and ambiguous codons at tips

- **IQ-TREE 3** makes the cell of a tip stop codon unknown and warns that it is "treated as missing data". A codon with any ambiguous or gap character becomes fully unknown, so `ACN` loses the information that it codes Thr ([alignment.cpp#L3830-L3870](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/alignment/alignment.cpp#L3830-L3870))
- **PAML** removes every column with a stop codon for all sequences, and expands `-` to all four bases ([treesub.c#L1004-L1040](https://github.com/abacus-gene/paml/blob/4c7902fe972737ef5e80bb18f159a6e6acace3d3/src/treesub.c#L1004-L1040))
- **HyPhy** stops when the alignment contains a stop codon, or, in its warning variant, deletes every column with a stop ([alignments.bf#L389-L422](https://github.com/veg/hyphy/blob/1fa8fd0469328b8d6df2c7273077170765beee5f/res/TemplateBatchFiles/libv3/tasks/alignments.bf#L389-L422))

All three erase the stop codons that tips share with their ancestors.

### Branch-length units of codon partitions

A codon model normalized to one substitution per codon needs branch lengths in substitutions per codon, which is 3 times the length per nucleotide site. IQ-TREE multiplies the rates of codon partitions by 3 ([phylosupertreeplen.cpp#L60-L61](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/tree/phylosupertreeplen.cpp#L60-L61)). HyPhy multiplies by 3 and fits a shared scale factor ([MG_REV.bf#L128](https://github.com/veg/hyphy/blob/1fa8fd0469328b8d6df2c7273077170765beee5f/res/TemplateBatchFiles/libv3/models/codon/MG_REV.bf#L128)).

## Translation rules

### Nextclade, Biopython and augur

Nextclade translates each codon on its own ([translate.rs#L8-L188](https://github.com/nextstrain/nextclade/blob/6190424cf8014d999fccf59d049c2c78a8bc28bb/packages/nextclade/src/translate/translate.rs#L8)):

- An ambiguous codon whose resolutions all code one amino acid gives that amino acid: `ACN` gives T, `YTR` gives L, `TAR` gives `*`
- Any other ambiguous codon, and any codon with one or two gaps, gives `X`
- `---` gives `-`, and a stop codon gives `*`

Biopython and augur differ in two ways:

- **Codons with several meanings**: Biopython writes `B`, `Z` or `J` (for example `RAY`, Asn or Asp, gives `B`), where Nextclade writes `X`. Auspice ignores only `X` in amino-acid mutations, so `B`, `Z` and `J` show as real mutations ([entropy.js#L323-L328](https://github.com/nextstrain/auspice/blob/db211cc6b44296d92529844502a95a3e68710449/src/util/entropy.js#L323-L328))
- **augur's fallback**: once any codon of a CDS has a partial gap, augur takes a fallback path that turns every stop codon of that CDS into `X` ([translate.py#L39-L96](https://github.com/nextstrain/augur/blob/0dc610704848c3a472f7715a1161277ef9b625a4/augur/translate.py#L39-L96)). The translation of one codon then depends on the other codons of the CDS

Both translators give `T` for `ACN`.

### Frameshifted sequences

Nextclade masks the codons between a frameshifting indel and its compensating indel with `X` ([translate_genes.rs#L160-L300](https://github.com/nextstrain/nextclade/blob/6190424cf8014d999fccf59d049c2c78a8bc28bb/packages/nextclade/src/translate/translate_genes.rs#L160)). A translation of reference-coordinate columns codon by codon instead reports reference-frame amino acids for a protein whose frame has changed.

### The `data/translate` script

`data/translate` made the bundled peptide files with Biopython's `forward_table.get(codon, "X")`, so it maps every ambiguous codon to `X`, `GCN` (Ala) included. Its peptides differ from the Nextclade rules at 464 tip residues of `rsv/a/2000` and 7,161 of `mpox/clade-ii/2000`. In every case the peptide file has `X` and the Nextclade rule gives an amino acid.

## Annotations and reading frames

### Overlapping reading frames

An inventory of 120 annotations (117 in the Nextclade datasets repository at [a231e7b](https://github.com/nextstrain/nextclade_data/tree/a231e7badf4747b9665708f7f662079b6f94264e), plus the 3 distinct annotations in `data/`) found:

- **Shifted-frame overlaps on the same strand in 46 annotations**: SARS-CoV-2 N and ORF9b over 294 nucleotides, influenza PB1 and PB1-F2, PA and PA-X, NS1 and NEP, M1 and M2, HIV-1, mpox, RSV M2-1 and M2-2, hMPV, mumps
- **Opposite-strand overlaps in 6 annotations**: mpox has 17 pairs of 2 to 42 nucleotides
- **Same-frame nested CDS in 10 annotations**: for example the PA-X prefix of PA, RT p66 and p51, prM and M

In `data/`, `rsv/a` has M2-1 and M2-2, `mpox/clade-ii` has 34 overlapping CDS rows, and `sc2` has 3 overlapping `gene` rows.

### Geometry of frameshifts, editing and splicing

The coordinates of the segments encode the special cases. A -1 frameshift or RNA editing shows as segments that overlap by 1, 2 or 4 nucleotides (SARS-CoV-2 nsp12 reads nucleotide 13468 twice; ebola GP repeats nucleotide 6923), a +1 frameshift as a gap of 1 nucleotide (influenza PA-X), and splicing as a long gap (influenza M2 and NEP). The `exception=` attribute only labels the case. RNA editing inserts a nucleotide that the template does not have, so one site read twice only approximates the biology.

### Prior art for overlapping frames

No published method reconstructs ancestral states across overlapping frames. Methods for overlapping genes estimate selection on sequence pairs or alignments without a tree <a id="cite-11"></a>[Sabath et al. 2008](https://doi.org/10.1371/journal.pone.0003996) [[11](#ref-11)]; <a id="cite-12"></a>[Wei and Zhang 2015](https://doi.org/10.1093/gbe/evu294) [[12](#ref-12)]; <a id="cite-13"></a>[Nelson et al. 2020](https://doi.org/10.1093/molbev/msaa087) [[13](#ref-13)]. Sabath et al. also show that independent models of the two frames count the shared data twice and bias the estimates. An exact model of interlocking frames is context-dependent and was fitted by MCMC on sequence pairs <a id="cite-2b"></a>[Pedersen and Jensen 2001](https://doi.org/10.1093/oxfordjournals.molbev.a003859) [[2](#ref-2)]. Phylogenetic hidden Markov models annotate viral overlaps but reconstruct no states <a id="cite-14"></a>[McCauley and Hein 2006](https://doi.org/10.1093/bioinformatics/btl092) [[14](#ref-14)]. A site-independent codon model applies only to nucleotides that belong to exactly one codon.

### Decoding across overlapping frames

When each site keeps one nucleotide value, an overlap region can be decoded jointly for all its frames: choose the nucleotide string of the region that maximizes the sum of the per-codon decoding scores of every frame. The codon score of each frame depends on three consecutive genome positions (a minus-strand codon is the reverse complement of such a window), so the objective is a chain of width-3 factors. Viterbi with the last two nucleotides as state is then exact at $O(64 \cdot L \cdot F)$ per node, for $L$ sites and $F$ frames.

- **Checked**: the dynamic program matched brute force in 200 of 200 random cases with shifted frames, opposite strands and three frames
- **Per-site argmax in overlap regions** is below the optimum in 194 of the same 200 cases; random posteriors overstate the real rate
- **The -1 frameshift fits the chain**: in `data/sc2`, ORF1a and ORF1b share nucleotide 13468 in two frames, and both codons are still windows
- **A spliced codon breaks the chain**: a codon split across a long gap (influenza M2, NEP) is not a window. It needs state carried across the junction or a small junction tree. No test covered this case

### Annotation reader of TreeTime

`fn parse_gff3_cds_features()` ([packages/treetime-io/src/gff.rs#L40-L136](../../packages/treetime-io/src/gff.rs#L40-L136)):

- **Reads `CDS` rows only**, so annotations with `gene` rows only give no CDS: SARS-CoV-2 `orfs`, dengue, yellow fever, and most influenza HA and NA files (38 of the 120). `data/sc2/*/annotation.gff3` has 12 `gene` rows and no `CDS` row
- **Groups rows by the first name attribute**, `Name` before `ID`, and ignores `Parent`. GFF3 joins the rows of one feature by a shared `ID`, and `Name` need not be unique ([GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/fe73505276dd324bf6a55773f3413fe2bed47af4/gff3.md?plain=1))
- **Rejects strand `.` and `?`**; every RSV row uses `.`
- **Ignores the phase column** ([kb/issues/M-ancestral-node-data-annotation-drops-cds-phase.md](../issues/M-ancestral-node-data-annotation-drops-cds-phase.md)), and sorts segments by coordinate. Its only checks are that each row has $1 \le \text{start} \le \text{end}$ and that the CDS length is a multiple of 3

Nextclade, which produces the annotations of the Nextstrain datasets, groups adjacent rows that share an `ID`, builds the hierarchy from `Parent`, uses a `gene` without `CDS` as a CDS, reads `.` and `?` as forward, computes the phase from the cumulative length, keeps segments in file order, and wraps circular genomes ([gene.rs#L42-L130](https://github.com/nextstrain/nextclade/blob/6190424cf8014d999fccf59d049c2c78a8bc28bb/packages/nextclade/src/gene/gene.rs#L42), [cds.rs#L46-L90](https://github.com/nextstrain/nextclade/blob/6190424cf8014d999fccf59d049c2c78a8bc28bb/packages/nextclade/src/gene/cds.rs#L46-L90)).

### Coordinate contract

GFF3 coordinates are reference coordinates, and nucleotide mutations are positions in alignment columns. The two agree only when the alignment columns are reference coordinates, with insertions relative to the reference removed. TreeTime compares the annotation with nothing in the alignment: the amino-acid route only checks that each peptide has one third of the CDS length (`fn validate_cds_length()` [packages/treetime/src/ancestral/aa.rs#L112-L127](../../packages/treetime/src/ancestral/aa.rs#L112-L127)), and the node data carries the annotation coordinates next to mutations in alignment columns. With an alignment that has insertions relative to the reference, the annotation and the mutations of the node data then use different coordinates, and a decoder that cuts codons from alignment columns would shift every codon after the first insertion, both without a message. The bundled annotated datasets fit: on `rsv/a/2000` and `mpox/clade-ii/2000` the alignment length equals the reference length, and every CDS is one row with phase 0. User data is not checked.

### Stop codons in the bundled data

The annotated terminal stop codon is not a sense state, and its position varies. In `rsv/a/20`, only 11 of 20 sequences of RSV-A G end at the annotated stop codon, 2 have an internal stop codon, and the others read past the annotated stop. No rule may assume a shared terminal stop codon.

### CDS names

The augur schema for root sequences allows CDS names that match `^[a-zA-Z0-9*_-]+$` ([schema-export-root-sequence.json#L25](https://github.com/nextstrain/augur/blob/0dc610704848c3a472f7715a1161277ef9b625a4/augur/data/schema-export-root-sequence.json#L25)), while the augur annotation schema also allows `.`, `(` and `)` ([schema-annotations.json#L31](https://github.com/nextstrain/augur/blob/0dc610704848c3a472f7715a1161277ef9b625a4/augur/data/schema-annotations.json#L31)). A name such as `ORF1a.1` passes one schema and fails the other. This is a defect of augur; a tool that writes both files can report such names.

## Evidence on codon substitution models

### Accuracy of ancestral states

- **The one direct benchmark is a tie.** On an experimental phylogeny with known ancestors, the per-site accuracy was 98.4% for the nucleotide model HKY+Γ and 98.3% for the codon model M0 with F3x4 frequencies <a id="cite-1b"></a>[Randall et al. 2016](https://doi.org/10.1038/ncomms12847) [[1](#ref-1)]
- **Codon models fit better**: GY94 fit 119 of 177 RNA-virus genes better than every nucleotide model <a id="cite-15"></a>[Shapiro et al. 2006](https://doi.org/10.1093/molbev/msj021) [[15](#ref-15)]
- **Nucleotide models recover recent divergences well.** The reported loss of 30% (deep) and 66% (recent) measures information lost in relative divergence-time estimates when an amino-acid model replaces a codon model <a id="cite-16"></a>[Ren et al. 2005](https://doi.org/10.1080/10635150500354688) [[16](#ref-16)]
- **Other work points the same way**: a study of ancestral protein estimation translates the per-site most probable DNA of a codon-position GTR model and compares no codon model <a id="cite-17"></a>[Hall 2006](https://doi.org/10.1073/pnas.0508991103) [[17](#ref-17)]; codon-partitioned nucleotide models reproduce the codon-model counts of nonsynonymous and synonymous changes within about 5% at 15 to 40 times less run time <a id="cite-18"></a>[Douglas and Bromham 2025](https://doi.org/10.64898/2025.12.21.695861) [[18](#ref-18)] (preprint); phylogenetic signal, not model heterogeneity, decides reconstruction accuracy <a id="cite-19"></a>[Muñiz-Trejo et al. 2025](https://doi.org/10.1093/molbev/msaf084) [[19](#ref-19)]; better-fitting models give reconstructions closer to held-out sequences at the protein level <a id="cite-20"></a>[Sennett and Theobald 2024](https://doi.org/10.1007/s00239-024-10162-3) [[20](#ref-20)]
- **An amino-acid process observed through codons is not Markovian** <a id="cite-21"></a>[Kosiol and Goldman 2011](https://doi.org/10.1016/j.jmb.2011.06.005) [[21](#ref-21)], which argues for modelling at the codon level rather than with empirical amino-acid matrices

> [!IMPORTANT]
> **Investigation required.** Three further accuracy studies compare no codon model in their abstracts; their full texts were not read <a id="cite-22"></a>[Aadland and Kolaczkowski 2020](https://doi.org/10.1093/gbe/evaa164) [[22](#ref-22)]; <a id="cite-23"></a>[Vialle et al. 2018](https://doi.org/10.1093/molbev/msy055) [[23](#ref-23)]; <a id="cite-24"></a>[Groussin et al. 2015](https://doi.org/10.1093/molbev/msu305) [[24](#ref-24)]. The unit of the 98.4% figure (Supplementary Table 2 of Randall et al.) was not read. FastML codon output, Bloom-lab tools, BEAST 2 codon packages and Lazarus were not examined.

### Model families

- **GY94 and MG94 in their original form** differ from the forms in common use. The original GY94 has no ω and weights changes by Grantham distance <a id="cite-25a"></a>[Goldman and Yang 1994](https://doi.org/10.1093/oxfordjournals.molbev.a040153) [[25](#ref-25)]; codeml uses the <a id="gloss-use-11"></a>M0 form <sup>[11](#gloss-11)</sup>. The original MG94 has separate synonymous and nonsynonymous rates and no κ <a id="cite-26a"></a>[Muse and Gaut 1994](https://doi.org/10.1093/oxfordjournals.molbev.a040152) [[26](#ref-26)]
- **<a id="gloss-use-12"></a>MG94xREV <sup>[12](#gloss-12)</sup>** uses GTR exchangeabilities, which describe skewed viral spectra. Expected SARS-CoV-2 counts at four-fold degenerate sites range from about 565 for C to T down to about 9 for T to G, which one transition-transversion ratio κ cannot describe <a id="cite-27a"></a>[Bloom and Neher 2023](https://doi.org/10.1093/ve/vead055) [[27](#ref-27)]. MG94xREV is the workhorse model of HyPhy ([MG_REV.bf#L22-L35](https://github.com/veg/hyphy/blob/1fa8fd0469328b8d6df2c7273077170765beee5f/res/TemplateBatchFiles/libv3/models/codon/MG_REV.bf#L22-L35))
- **MG94xREV nests the codon-position nucleotide model**: with ω = 1 and a genetic code without stop codons, its generator is the <a id="gloss-use-13"></a>Kronecker sum <sup>[13](#gloss-13)</sup> of three per-position GTR generators with shared exchangeabilities, so the two models must give equal likelihoods
- **Bayesian model comparison** favored MG models with codon or amino-acid preference terms and Dirichlet-process rate variation; plain MG94 was not compared with plain GY94, and GY-type models with F61 fit better on the HIV set. In GY94 the rate of a change depends on the two unchanged nucleotides, which the authors call undesirable <a id="cite-28"></a>[Rodrigue et al. 2008](https://doi.org/10.1534/genetics.108.092254) [[28](#ref-28)]

### Codon frequencies

- **F3x4 biases the frequencies** that the model induces once stop codons are removed, and the GTR rate estimates. <a id="gloss-use-14"></a>CF3x4 <sup>[14](#gloss-14)</sup> corrects this with a solver for 9 unknowns <a id="cite-29a"></a>[Kosakovsky Pond et al. 2010](https://doi.org/10.1371/journal.pone.0011230) [[29](#ref-29)]
- **MG94 with F61 frequencies** is not a separate model: MG94 uses per-position nucleotide frequencies, and the per-position marginals of the F61 frequencies are the F3x4 counts. IQ-TREE rejects the combination ([modelcodon.cpp#L550](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/modelcodon.cpp#L550))
- **Unobserved codons under F61** need a floor: IQ-TREE raises frequencies below 1e-4 ([alignment.cpp#L5882-L5904](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/alignment/alignment.cpp#L5882-L5904))

### Rate variation and multi-nucleotide changes

- **ω varies strongly among sites** in every data set studied, influenza HA and HIV-1 included <a id="cite-30"></a>[Yang et al. 2000](https://doi.org/10.1093/genetics/155.1.431) [[30](#ref-30)], so one global ω is a first model, not an adequate one
- **Rate categories among codons and instantaneous multi-nucleotide changes improve fit** <a id="cite-31"></a>[Whelan and Goldman 2004](https://doi.org/10.1534/genetics.103.023226) [[31](#ref-31)]; <a id="cite-32"></a>[Lucaci et al. 2021](https://doi.org/10.1371/journal.pone.0248337) [[32](#ref-32)]. Some multi-nucleotide mutations recur hundreds of times in SARS-CoV-2 <a id="cite-33"></a>[De Maio et al. 2025](https://doi.org/10.1093/molbev/msaf272) [[33](#ref-33)]. A single-change codon model needs two or three separate events for such a change and penalizes it; a nucleotide model treats the changes as independent. No source shows an effect on ancestral states
- **Viral mutation spectra are not reversible**: the SARS-CoV-2 spectrum <a id="cite-27b"></a>[Bloom and Neher 2023](https://doi.org/10.1093/ve/vead055) [[27](#ref-27)] cannot be captured by a reversible codon or nucleotide model. The effect on ancestral states is unknown

> [!IMPORTANT]
> **Investigation required.** How well ω is identified on near-clonal data such as `data/*/20` is not quantified; wide confidence intervals are expected.

## Cost and representation

- **Decoding is cheap**: 64 products and 64 comparisons per codon per node, on site posteriors that the dense and the sparse representation already hold. It was not measured
- **Dense 61-state storage does not fit.** The dense path keeps four sites-by-states arrays of 64-bit floats per edge after the passes ([packages/treetime/src/partition/storage/dense.rs#L11-L78](../../packages/treetime/src/partition/storage/dense.rs#L11-L78)). At 61 states this is about 26 GB for `dengue/2000`, 36 GB for `rsv/a/2000`, 195 GB for `sc2/5100` and 430 GB for `mpox/clade-ii/2000` (estimates from the array shapes). Sparse is the default representation ([packages/treetime/src/partition/algo/infer_dense.rs#L1-L7](../../packages/treetime/src/partition/algo/infer_dense.rs#L1-L7)). A dense codon analysis of 997 dengue genomes needed 21.5 GB on GPUs <a id="cite-34"></a>[Ayres et al. 2019](https://doi.org/10.1093/sysbio/syz020) [[34](#ref-34)]
- **Sparse 61-state runs are estimated at 0.2 to 4 GB** from array shapes and the variation within subtrees (the median internal node has 8 to 53 variable codon columns in its subtree). The sparse nucleotide run of `mpox/clade-ii/2000` alone peaks at 2.59 GB, so the estimate of 1.1 GB for an mpox codon partition is doubtful. An extrapolation by the square of the state count from the 21-state route gives about 6 single-thread minutes for a 61-state mpox run
- **Identical codon columns are rare**: 9,787 of the 9,814 codon columns of `sc2/5100` are unique, so collapsing duplicates saves almost nothing
- **Three sparse costs grow with the state count**: the shared per-state profiles of each message (about 33 KB per message, about 1.7 GB for `sc2/5100`), one heap allocation per stored site, and the transition counting, which allocates an array of L by n by n per edge (`fn get_branch_mutation_matrix()` [packages/treetime/src/gtr/infer_gtr/common.rs#L132-L194](../../packages/treetime/src/gtr/infer_gtr/common.rs#L132-L194); 0.29 GB per edge for `sc2/5100` and 1.65 GB for `mpox/clade-ii/2000` at 61 states). The sparse threshold for resolved sites is fixed at 1e-4 ([packages/treetime/src/partition/marginal/sparse/message.rs#L20](../../packages/treetime/src/partition/marginal/sparse/message.rs#L20)) and untested at 61 states
- **The web build is capped at 4 GiB** (wasm32, `MAXIMUM_MEMORY=4GB`, [.cargo/config.toml#L74-L75](../../.cargo/config.toml#L74-L75)). Sparse codon runs up to about 2,000 tips should fit; `sc2/5100` is at the limit. A memory64 build allows 16 GB, but it does not run in Safari and is 10% to 2 times slower

> [!IMPORTANT]
> **Investigation required.** The sparse 61-state memory and time figures are estimates. Measure them on a codon partition before any design depends on them.

## Sampling

`--sample-from-profile` draws the state of each node from that node's own marginal posterior (`fn sample_from_profile()` [packages/treetime/src/partition/marginal/sample.rs#L39-L55](../../packages/treetime/src/partition/marginal/sample.rs#L39-L55), called with each node's profile in the dense and the sparse partition). The draws of a parent and its child are independent, so a sampled set of node states is not a sample of a history, and the mutations between sampled nodes are not sampled mutations. This also affects nucleotide runs without an annotation. A draw from the joint posterior needs pre-order sampling, as BEAST does.

## Pitfalls found in an unmerged prototype

An earlier prototype added GY94 and MG94 codon models, a codon partition, and an overlay that wrote codon results over a nucleotide reconstruction of the whole alignment. It is not merged, and it is written against graph and partition interfaces that no longer exist. Its defects show where a codon implementation goes wrong:

- **Branches 3 times too short**: the codon matrix was normalized to one substitution per codon but received branch lengths per nucleotide site (section [Branch-length units of codon partitions](#branch-length-units-of-codon-partitions)). Histories with $k$ extra changes then lose weight by about $3^k$, so reported confidences are too high, and the most probable codon can change where histories trade the number of changes against their kind (derived, not measured)
- **`MG94 + F61` equal to MG94 + F3x4**: the option read the per-position marginals of F61, which are the F3x4 counts
- **Overlapping CDS rejected**: this fails on every annotated dataset in `data/` and on 46 of the 120 virus annotations
- **Leaf rows of 1/k**: an ambiguous tip codon got 1/k on each of its k compatible codons, while the nucleotide path uses indicator rows (1.0). Rows of 1/k shift the log-likelihood by $\sum \log(1/k)$ per leaf codon, so codon and nucleotide likelihoods are not comparable
- **Tip stop codons erased**: tip stop codons became missing data, so the 61-state model wrote sense codons at internal nodes whose tips share a stop codon
- **Split decoder**: nucleotides came from the most probable codon, amino acids from the summed amino-acid posterior. Nucleotide and amino-acid outputs disagreed at 0.17% of residues on `rsv/a/20`, and the prototype's own test asserted codon F against amino acid L
- **Coding data in two models**: the nucleotide reconstruction of the whole alignment and the codon overlay both used the coding sites
- **Codon model names in the global model enum**: commands without codon support accepted the names and failed at run time
- **A self-confirming oracle**: the "analytical oracle" of the matrix tests reused the code under test

The standard code table, the GY94 and MG94 matrix forms, and the triplet expansion with reverse complement were correct.

## Claims that the evidence refutes

| Claim | Evidence |
| --- | --- |
| In codon space the inconsistency of nucleotide and amino-acid outputs cannot arise | A codon model with a split decoder gives inconsistent outputs (section [Decoding rules from a loss function](#decoding-rules-from-a-loss-function)); consistency comes from translating one decoded codon |
| The order of reconstruction and translation (argmax then translation, against translation then reconstruction) is the main source of disagreement between the routes | In v1, model differences at confident codons and gap handling cause most differences (section [Results](#results)) |
| augur uses JC69 for nucleotides and a uniform 22-state JC69 for amino acids | `augur ancestral` infers a GTR model after a JC69 seed, for both alphabets |
| For `ACN`, Nextclade writes `T` and augur writes `X` | Both write `T`; the translators differ at codons with several meanings and in augur's fallback path |
| The RSV pipeline uses the peptide route | It runs both routes, and the translation route wins the export merge |
| No published benchmark compares codon and nucleotide ancestral states | Randall et al. 2016 compare them: a tie |
| The 30% and 66% losses of Ren et al. 2005 measure phylogenetic information for topology | They measure information lost in relative divergence-time estimates |
| Tip stop codons are uninformative, and mapping them to missing data is harmless | Real stop codons exist at internal nodes |
| F3x4 is a robust default for codon frequencies | F3x4 biases induced frequencies and rate estimates; CF3x4 corrects this |

## Glossary

1. <a id="gloss-1"></a> **Diffuse codon.** A codon at one node whose most probable triplet has a posterior below 1/2. Only at such codons can the decoding rules choose different codons. [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Sense codon.** One of the 61 triplets that code an amino acid in the standard genetic code. The other 3 triplets (TAA, TAG, TGA) are stop codons. [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **Per-site argmax.** The most probable nucleotide at each alignment column of a node, chosen column by column. TreeTime writes this nucleotide in its marginal reconstruction today. [↩](#gloss-use-3)
4. <a id="gloss-4"></a> **Bayes rule.** The decision that minimizes the posterior expected loss. For 0-1 loss it is the most probable state. [↩](#gloss-use-4)
5. <a id="gloss-5"></a> **Lexicographic loss.** A loss that ranks one error kind strictly above another: any amino-acid error costs more than every nucleotide error. It is the limit of the weighted loss for a large amino-acid weight. [↩](#gloss-use-5)
6. <a id="gloss-6"></a> **Centroid estimator.** The estimate that minimizes the expected Hamming loss. For independent components it takes the most probable value of each component <a id="cite-3b"></a>[Carvalho and Lawrence 2008](https://doi.org/10.1073/pnas.0712329105) [[3](#ref-3)]. [↩](#gloss-use-6)
7. <a id="gloss-7"></a> **Joint reconstruction.** The single most probable assignment of states to all internal nodes together, computed by a dynamic program over the tree <a id="cite-4b"></a>[Pupko et al. 2000](https://doi.org/10.1093/oxfordjournals.molbev.a026369) [[4](#ref-4)]. Marginal reconstruction instead gives each node its own posterior over states, from all tip data <a id="cite-35"></a>[Yang et al. 1995](https://doi.org/10.1093/genetics/141.4.1641) [[35](#ref-35)]; <a id="cite-36"></a>[Koshi and Goldstein 1996](https://doi.org/10.1007/bf02198858) [[36](#ref-36)]. [↩](#gloss-use-7)
8. <a id="gloss-8"></a> **Marginal MAP.** Maximization over some variables of a sum over the others, here the most probable amino-acid history summed over synonymous codons. It is NP-hard for general graphical models <a id="cite-37"></a>[Park and Darwiche 2004](https://doi.org/10.1613/jair.1236) [[37](#ref-37)]. [↩](#gloss-use-8)
9. <a id="gloss-9"></a> **Robust counting.** Expected numbers of labelled substitutions (for example synonymous and nonsynonymous) on each branch, computed over the posterior of histories instead of from one reconstruction <a id="cite-8b"></a>[Minin and Suchard 2008](https://doi.org/10.1007/s00285-007-0120-8) [[8](#ref-8)]; <a id="cite-9b"></a>[Lemey et al. 2012](https://doi.org/10.1093/bioinformatics/bts580) [[9](#ref-9)]. [↩](#gloss-use-9)
10. <a id="gloss-10"></a> **Pre-order sampling.** A draw from the joint posterior of all ancestral states: the root from its posterior, then each child from its conditional distribution given the drawn parent state and the data below the child. [↩](#gloss-use-10)
11. <a id="gloss-11"></a> **M0 form of GY94.** The simplified Goldman and Yang codon model that codeml uses: rate $\pi_j$, times κ for transitions and ω for nonsynonymous changes. The original model weights changes by Grantham distance and has no ω <a id="cite-25b"></a>[Goldman and Yang 1994](https://doi.org/10.1093/oxfordjournals.molbev.a040153) [[25](#ref-25)]. [↩](#gloss-use-11)
12. <a id="gloss-12"></a> **MG94xREV.** The Muse and Gaut codon model <a id="cite-26b"></a>[Muse and Gaut 1994](https://doi.org/10.1093/oxfordjournals.molbev.a040152) [[26](#ref-26)] with GTR nucleotide exchangeabilities, a nonsynonymous rate ratio ω, and the frequency of the target nucleotide in each rate. [↩](#gloss-use-12)
13. <a id="gloss-13"></a> **Kronecker sum.** The generator of three independent processes on one combined state space, $Q_1 \oplus Q_2 \oplus Q_3$. Its transition matrix is the Kronecker product of the three transition matrices. [↩](#gloss-use-13)
14. <a id="gloss-14"></a> **F3x4 and CF3x4.** Codon frequencies as products of nucleotide frequencies per codon position. F3x4 uses the observed frequencies directly, which biases the frequencies that the model induces once stop codons are removed. CF3x4 solves for corrected nucleotide frequencies that reproduce the observed ones <a id="cite-29b"></a>[Kosakovsky Pond et al. 2010](https://doi.org/10.1371/journal.pone.0011230) [[29](#ref-29)]. [↩](#gloss-use-14)

## References

1. <a id="ref-1"></a> Randall, Ryan N., Caelan E. Radford, Kelsey A. Roof, Divya K. Natarajan, and Eric A. Gaucher. 2016. "An experimental phylogeny to benchmark ancestral sequence reconstruction." _Nature Communications_ 7:12847. https://doi.org/10.1038/ncomms12847 [↩¹](#cite-1a) [↩²](#cite-1b)
2. <a id="ref-2"></a> Pedersen, Anne-Mette Krabbe, and Jens Ledet Jensen. 2001. "A dependent-rates model and an MCMC-based methodology for the maximum-likelihood analysis of sequences with overlapping reading frames." _Molecular Biology and Evolution_ 18:763-776. https://doi.org/10.1093/oxfordjournals.molbev.a003859 [↩¹](#cite-2a) [↩²](#cite-2b)
3. <a id="ref-3"></a> Carvalho, Luis E., and Charles E. Lawrence. 2008. "Centroid estimation in discrete high-dimensional spaces with applications in biology." _Proceedings of the National Academy of Sciences_ 105:3209-3214. https://doi.org/10.1073/pnas.0712329105 [↩¹](#cite-3a) [↩²](#cite-3b)
4. <a id="ref-4"></a> Pupko, Tal, Itsik Pe'er, Ron Shamir, and Dan Graur. 2000. "A fast algorithm for joint reconstruction of ancestral amino acid sequences." _Molecular Biology and Evolution_ 17:890-896. https://doi.org/10.1093/oxfordjournals.molbev.a026369 [↩¹](#cite-4a) [↩²](#cite-4b)
5. <a id="ref-5"></a> Williams, Paul D., David D. Pollock, Benjamin P. Blackburne, and Richard A. Goldstein. 2006. "Assessing the accuracy of ancestral protein reconstruction methods." _PLoS Computational Biology_ 2:e69. https://doi.org/10.1371/journal.pcbi.0020069 [↩](#cite-5)
6. <a id="ref-6"></a> Krishnan, Neeraja M., Hervé Seligmann, Caro-Beth Stewart, A. P. Jason de Koning, and David D. Pollock. 2004. "Ancestral sequence reconstruction in primate mitochondrial DNA: compositional bias and effect on functional inference." _Molecular Biology and Evolution_ 21:1871-1883. https://doi.org/10.1093/molbev/msh198 [↩](#cite-6)
7. <a id="ref-7"></a> Matsumoto, Tomotaka, Hiroshi Akashi, and Ziheng Yang. 2015. "Evaluation of ancestral sequence reconstruction methods to infer nonstationary patterns of nucleotide substitution." _Genetics_ 200:873-890. https://doi.org/10.1534/genetics.115.177386 [↩](#cite-7)
8. <a id="ref-8"></a> Minin, Vladimir N., and Marc A. Suchard. 2008. "Counting labeled transitions in continuous-time Markov models of evolution." _Journal of Mathematical Biology_ 56:391-412. https://doi.org/10.1007/s00285-007-0120-8 [↩¹](#cite-8a) [↩²](#cite-8b)
9. <a id="ref-9"></a> Lemey, Philippe, Vladimir N. Minin, Filip Bielejec, Sergei L. Kosakovsky Pond, and Marc A. Suchard. 2012. "A counting renaissance: combining stochastic mapping and empirical Bayes to quickly detect amino acid sites under positive selection." _Bioinformatics_ 28:3248-3256. https://doi.org/10.1093/bioinformatics/bts580 [↩¹](#cite-9a) [↩²](#cite-9b)
10. <a id="ref-10"></a> Oliva, A., S. Pulicani, V. Lefort, L. Bréhélin, O. Gascuel, and S. Guindon. 2019. "Accounting for ambiguity in ancestral sequence reconstruction." _Bioinformatics_ 35:4290-4297. Preprint: https://doi.org/10.1101/409029 https://doi.org/10.1093/bioinformatics/btz249 [↩](#cite-10)
11. <a id="ref-11"></a> Sabath, Niv, Giddy Landan, and Dan Graur. 2008. "A method for the simultaneous estimation of selection intensities in overlapping genes." _PLoS ONE_ 3:e3996. https://doi.org/10.1371/journal.pone.0003996 [↩](#cite-11)
12. <a id="ref-12"></a> Wei, Xinzhu, and Jianzhi Zhang. 2015. "A simple method for estimating the strength of natural selection on overlapping genes." _Genome Biology and Evolution_ 7:381-390. https://doi.org/10.1093/gbe/evu294 [↩](#cite-12)
13. <a id="ref-13"></a> Nelson, Chase W., Zachary Ardern, and Xinzhu Wei. 2020. "OLGenie: estimating natural selection to predict functional overlapping genes." _Molecular Biology and Evolution_ 37:2440-2449. https://doi.org/10.1093/molbev/msaa087 [↩](#cite-13)
14. <a id="ref-14"></a> McCauley, Stephen, and Jotun Hein. 2006. "Using hidden Markov models and observed evolution to annotate viral genomes." _Bioinformatics_ 22:1308-1316. https://doi.org/10.1093/bioinformatics/btl092 [↩](#cite-14)
15. <a id="ref-15"></a> Shapiro, Beth, Andrew Rambaut, and Alexei J. Drummond. 2006. "Choosing appropriate substitution models for the phylogenetic analysis of protein-coding sequences." _Molecular Biology and Evolution_ 23:7-9. https://doi.org/10.1093/molbev/msj021 [↩](#cite-15)
16. <a id="ref-16"></a> Ren, Fengrong, Hiroshi Tanaka, and Ziheng Yang. 2005. "An empirical examination of the utility of codon-substitution models in phylogeny reconstruction." _Systematic Biology_ 54:808-818. https://doi.org/10.1080/10635150500354688 [↩](#cite-16)
17. <a id="ref-17"></a> Hall, Barry G. 2006. "Simple and accurate estimation of ancestral protein sequences." _Proceedings of the National Academy of Sciences_ 103:5431-5436. https://doi.org/10.1073/pnas.0508991103 [↩](#cite-17)
18. <a id="ref-18"></a> Douglas, Jordan, and Lindell Bromham. 2025. "Reconstructing substitution histories on phylogenies, with accuracy, precision, and coverage." _bioRxiv_ preprint. https://doi.org/10.64898/2025.12.21.695861 [↩](#cite-18)
19. <a id="ref-19"></a> Muñiz-Trejo, Ricardo, Yeonwoo Park, and Joseph W. Thornton. 2025. "Robustness of ancestral sequence reconstruction to among-site and among-lineage evolutionary heterogeneity." _Molecular Biology and Evolution_ 42:msaf084. https://doi.org/10.1093/molbev/msaf084 [↩](#cite-19)
20. <a id="ref-20"></a> Sennett, Michael A., and Douglas L. Theobald. 2024. "Extant sequence reconstruction: the accuracy of ancestral sequence reconstructions evaluated by extant sequence cross-validation." _Journal of Molecular Evolution_ 92:181-206. https://doi.org/10.1007/s00239-024-10162-3 [↩](#cite-20)
21. <a id="ref-21"></a> Kosiol, Carolin, and Nick Goldman. 2011. "Markovian and non-Markovian protein sequence evolution: aggregated Markov process models." _Journal of Molecular Biology_ 411:910-923. https://doi.org/10.1016/j.jmb.2011.06.005 [↩](#cite-21)
22. <a id="ref-22"></a> Aadland, Kelsey, and Bryan Kolaczkowski. 2020. "Alignment-integrated reconstruction of ancestral sequences improves accuracy." _Genome Biology and Evolution_ 12:1549-1565. https://doi.org/10.1093/gbe/evaa164 [↩](#cite-22)
23. <a id="ref-23"></a> Vialle, Ricardo Assunção, Asif U. Tamuri, and Nick Goldman. 2018. "Alignment modulates ancestral sequence reconstruction accuracy." _Molecular Biology and Evolution_ 35:1783-1797. https://doi.org/10.1093/molbev/msy055 [↩](#cite-23)
24. <a id="ref-24"></a> Groussin, Mathieu, Joanne K. Hobbs, Gergely J. Szöllősi, Simonetta Gribaldo, Vickery L. Arcus, and Manolo Gouy. 2015. "Toward more accurate ancestral protein genotype-phenotype reconstructions with the use of species tree-aware gene trees." _Molecular Biology and Evolution_ 32:13-22. https://doi.org/10.1093/molbev/msu305 [↩](#cite-24)
25. <a id="ref-25"></a> Goldman, Nick, and Ziheng Yang. 1994. "A codon-based model of nucleotide substitution for protein-coding DNA sequences." _Molecular Biology and Evolution_ 11:725-736. https://doi.org/10.1093/oxfordjournals.molbev.a040153 [↩¹](#cite-25a) [↩²](#cite-25b)
26. <a id="ref-26"></a> Muse, Spencer V., and Brandon S. Gaut. 1994. "A likelihood approach for comparing synonymous and nonsynonymous nucleotide substitution rates, with application to the chloroplast genome." _Molecular Biology and Evolution_ 11:715-724. https://doi.org/10.1093/oxfordjournals.molbev.a040152 [↩¹](#cite-26a) [↩²](#cite-26b)
27. <a id="ref-27"></a> Bloom, Jesse D., and Richard A. Neher. 2023. "Fitness effects of mutations to SARS-CoV-2 proteins." _Virus Evolution_ 9:vead055. https://doi.org/10.1093/ve/vead055 [↩¹](#cite-27a) [↩²](#cite-27b)
28. <a id="ref-28"></a> Rodrigue, Nicolas, Nicolas Lartillot, and Hervé Philippe. 2008. "Bayesian comparisons of codon substitution models." _Genetics_ 180:1579-1591. https://doi.org/10.1534/genetics.108.092254 [↩](#cite-28)
29. <a id="ref-29"></a> Kosakovsky Pond, Sergei, Wayne Delport, Spencer V. Muse, and Konrad Scheffler. 2010. "Correcting the bias of empirical frequency parameter estimators in codon models." _PLoS ONE_ 5:e11230. https://doi.org/10.1371/journal.pone.0011230 [↩¹](#cite-29a) [↩²](#cite-29b)
30. <a id="ref-30"></a> Yang, Ziheng, Rasmus Nielsen, Nick Goldman, and Anne-Mette Krabbe Pedersen. 2000. "Codon-substitution models for heterogeneous selection pressure at amino acid sites." _Genetics_ 155:431-449. https://doi.org/10.1093/genetics/155.1.431 [↩](#cite-30)
31. <a id="ref-31"></a> Whelan, Simon, and Nick Goldman. 2004. "Estimating the frequency of events that cause multiple-nucleotide changes." _Genetics_ 167:2027-2043. https://doi.org/10.1534/genetics.103.023226 [↩](#cite-31)
32. <a id="ref-32"></a> Lucaci, Alexander G., Sadie R. Wisotsky, Stephen D. Shank, Steven Weaver, and Sergei L. Kosakovsky Pond. 2021. "Extra base hits: widespread empirical support for instantaneous multiple-nucleotide changes." _PLoS ONE_ 16:e0248337. https://doi.org/10.1371/journal.pone.0248337 [↩](#cite-32)
33. <a id="ref-33"></a> De Maio, Nicola, Olivier Anoufa, Kyle Smith, Yatish Turakhia, and Nick Goldman. 2025. "Highly recurrent multinucleotide mutations in SARS-CoV-2." _Molecular Biology and Evolution_ 42:msaf272. https://doi.org/10.1093/molbev/msaf272 [↩](#cite-33)
34. <a id="ref-34"></a> Ayres, Daniel L., Michael P. Cummings, Guy Baele, et al. 2019. "BEAGLE 3: improved performance, scaling, and usability for a high-performance computing library for statistical phylogenetics." _Systematic Biology_ 68:1052-1061. https://doi.org/10.1093/sysbio/syz020 [↩](#cite-34)
35. <a id="ref-35"></a> Yang, Ziheng, Sudhir Kumar, and Masatoshi Nei. 1995. "A new method of inference of ancestral nucleotide and amino acid sequences." _Genetics_ 141:1641-1650. https://doi.org/10.1093/genetics/141.4.1641 [↩](#cite-35)
36. <a id="ref-36"></a> Koshi, Jeffrey M., and Richard A. Goldstein. 1996. "Probabilistic reconstruction of ancestral protein sequences." _Journal of Molecular Evolution_ 42:313-320. https://doi.org/10.1007/bf02198858 [↩](#cite-36)
37. <a id="ref-37"></a> Park, James D., and Adnan Darwiche. 2004. "Complexity results and approximation strategies for MAP explanations." _Journal of Artificial Intelligence Research_ 21:101-133. https://doi.org/10.1613/jair.1236 [↩](#cite-37)
