# Partition Configuration for TreeTime v1

## 1. Current state

Ancestral, optimize and timetree build their sequence partitions from partition plans. A plan is a set of sites (columns) of one `--alignment` file with a model description (`struct PartitionPlan`, [packages/treetime/src/partition/resolve.rs#L124-L127](../../packages/treetime/src/partition/resolve.rs#L124-L127)). Without `--partition` entries, each file is one plan with the model of `--model` and `--model-params`; with entries, each entry is one plan with its own sites and model. One resolver (`fn resolve_partitions()`) turns the entries and the alignment lengths into the plans for the CLI and the apps, and the pipelines split each file into the sequences of its partitions (`fn split_alignments()`, [packages/treetime/src/partition/split.rs#L57-L101](../../packages/treetime/src/partition/split.rs#L57-L101)). Prune builds one sparse JC69 partition per file from the layouts of the same resolver, and homoplasy one partition from its one file. [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md) records the rules.

Optimize and timetree hold their partitions in `enum PartitionSet` ([packages/treetime/src/partition/set.rs#L70](../../packages/treetime/src/partition/set.rs#L70)), all dense or all sparse as `--representation` selects; the default comes from `infer_dense()`, which always returns `false` (sparse). Results come back per alignment file (`struct AlignmentResult`) and per partition (`struct PartitionResult`) ([packages/treetime/src/partition/result.rs](../../packages/treetime/src/partition/result.rs)), and the partitions of a split file are assembled back to the coordinates of the file ([kb/algo/partition-assembly.md](../algo/partition-assembly.md)). Mugration runs one discrete trait partition alone.

v0 has no explicit partition concept. `TreeAnc` owns a single `self.gtr` and single `self.data`. One alignment, one model.

## 2. Configuration methods

Five distinct methods for specifying what partitions to create have been discussed in the project. They are not mutually exclusive: a partition configuration layer would accept input from any of these sources and produce a uniform `Vec<PartitionSpec>`.

### 2.1 `--representation` flag

The `--representation dense|sparse` flag of ancestral, optimize and timetree selects the representation of every partition of the run. Without the flag, `infer_dense()` decides (S2.2). Prune always uses the sparse representation with JC69.

**Sources:** `packages/app-commands/src/commands/shared/representation.rs`

### 2.2 Auto-detection heuristic (`infer_dense`)

When `--representation` is omitted, `infer_dense()` in `packages/treetime/src/partition/algo/infer_dense.rs` determines dense vs sparse. The stub always returns `false`. The planned heuristic: select dense when tree branches are long enough that most positions become variable (expected mutations per branch approaches sequence length). Below that threshold, sparse saves memory and compute.

**Sources:** `kb/issues/N-representation-infer-dense-stub.md`

### 2.3 Multiple alignment files and partition entries

Each `--alignment` file is numbered by its position in the list. Without `--partition` entries, each file is one partition and all files take the one `--model` choice. With entries, each entry takes a set of sites of one file (site ranges `a`, `a-b`, `a-b\k`, numbered from 1) and carries its own model, so one file can hold codon positions or genes with separate models. Optimize and timetree estimate one relative rate per partition. Discrete characters next to sequences still need a config format (S2.5).

**Sources:** [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md), `packages/app-commands/src/commands/shared/partition.rs`

### 2.4 Segmented genomes

For segmented viruses (influenza: 8 segments, same alphabet, independent evolution), segments that share one tree run as one `--alignment` file per segment, or as `--partition` entries of one file (S2.3).

Segments with different trees use `--segment-trees <tree>... --segment-alignments <fasta>...` of `ancestral`, `optimize` and `timetree`: one tree and one alignment per segment, merged by node name into one graph whose edges record which segments they carry. The bundled `data/flu/h3n2/segments/*/pipeline-nextstrain.yaml` run on such a graph. Open gaps are in [kb/issues/N-io-multi-segment-genome-input.md](../issues/N-io-multi-segment-genome-input.md) and [kb/issues/M-partition-passes-crossover-scopes-and-discrete-traits.md](../issues/M-partition-passes-crossover-scopes-and-discrete-traits.md).

v0 handles segments via concatenation.

**Sources:** `kb/issues/N-io-multi-segment-genome-input.md`, `kb/proposals/edge-scoped-inheritance.md`

### 2.5 Config file (YAML/JSON)

A full proposal exists for a structured config file specifying:

- Partition list: each entry has alignment path, alphabet, model name or model file, dense/sparse mode
- Shared parameters: tree path, clock model, coalescent settings
- Output specification

This handles all use cases: multi-segment genomes, mixed data types (nucleotide + discrete location), per-partition model assignment, model file references for two-stage estimation. The config file complements CLI flags (simple cases use CLI, complex cases use config). Per-partition models and site ranges of sequence partitions are implemented as the `partition` setting (S2.3), which a per-command config file can hold.

Prior art surveyed in the proposal:

| Tool     | Format                                   | Partition specification                         |
| -------- | ---------------------------------------- | ----------------------------------------------- |
| NEXUS    | SETS block                               | Character ranges within concatenated alignment  |
| IQ-TREE  | `-spp`/`-q`/`-sp` flags + partition file | Column ranges, three edge-linking modes         |
| RAxML-NG | Partition file                           | `MODEL, name = start-end` per line              |
| BEAST    | XML                                      | Full model hierarchy, link/unlink per partition |

Open design questions: config files that span several commands; edge-linked vs edge-unlinked partitions.

**Sources:** `kb/proposals/config-file-multi-partition.md`

### 2.6 Unified input formats

Auspice JSON and UShER MAT files contain embedded sequence data. Loading these formats could populate partitions directly without separate alignment files.

**Sources:** `kb/proposals/unified-input-format-support.md`

## 3. Automatic partitioning (`--partition-by`)

Partitions are derived from an annotation or from the alignment itself with `--partition-by` (`enum PartitionScheme`, [packages/app-commands/src/commands/shared/partition_by.rs#L28](../../packages/app-commands/src/commands/shared/partition_by.rs#L28)). A scheme produces ordinary partition entries with explicit site ranges, which the resolver of the entries then checks with its own rules and messages, so a scheme is a generator of entries and not a separate code path. `--partition-by` and `--partition` conflict. [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md) holds the full rules, including the annotation contract and the overlap ordering.

### 3.1 Rate-based partitioning (`rate-categories`)

A Fitch reconstruction on the input tree counts, for each site, the edges with a substitution at that site (`fn fitch_site_change_counts()`, [packages/treetime/src/partition/fitch/site_changes.rs#L11](../../packages/treetime/src/partition/fitch/site_changes.rs#L11)), and `--rate-category-rule` sorts the sites into `slow`, `medium` and `fast`: `change-counts` (default; 0, 1, 2 or more changes), `median-split`, or `k-means` (an exact one-dimensional k-means over the counts). The scheme needs no annotation.

The number of buckets is fixed at three. `rate_category_buckets()` maps the thresholds onto the name array `["slow", "medium", "fast"]` and `k_means_thresholds()` returns at most two split points ([packages/app-commands/src/commands/shared/partition_by.rs#L257-L292](../../packages/app-commands/src/commands/shared/partition_by.rs#L257-L292)), so there is no setting for the number of categories. A bucket that stays empty is left out. Among-site rate variation with a chosen number of categories is `--gamma-categories K`, a different model: it averages each site over K rates instead of assigning it to one bucket ([kb/decisions/discrete-gamma-rate-variation.md](../decisions/discrete-gamma-rate-variation.md)).

The buckets are recorded as biased in [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md), with a simulation measurement, and PartitionFinder 2 withdrew k-means splits of sites for DNA and protein data for the same reason. S3.3 measures the circularity.

### 3.2 Annotation-based partitioning (`codon-positions` and `genes`)

`codon-positions` puts codon positions 1, 2 and 3 of every CDS of `--annotation` into `cp1`, `cp2` and `cp3`, and the remaining sites into `noncoding`. `genes` makes one partition per CDS, named by the CDS, plus `noncoding`. Both read GFF3 through `fn gff3_cds_read()` ([packages/treetime-io/src/gff.rs](../../packages/treetime-io/src/gff.rs), on `noodles-gff`), take one annotation file per `--alignment` file, give an overlapped site to the longest CDS that covers it, and take the codon position of a site from the phase of its own segment, on both strands.

Both are implemented (`kb/features/representation.md` lists the schemes `[x]`). The annotation contract is strict: the CDS coordinates must be alignment columns, so each file needs a `##sequence-region` line that runs from 1 to the length of its alignment.

### 3.3 Measured behavior and recommended defaults

Measured with `treetime optimize` on `data/dengue/2000` (1954 genomes, 10723 sites, DENV2 reference coordinates), using the annotation of `data/dengue/20`. Baseline without partition settings: tree length 5.364.

| Scheme | Partitions | Relative rates mu | Tree length |
| --- | --- | --- | --- |
| `codon-positions` (default sharing: none) | cp1, cp2, cp3, noncoding | 0.421, 0.184, 2.382, 1.075 | 5.522 |
| `genes` (default sharing: W and pi) | 11 CDS plus noncoding | 0.843 to 1.128 | 5.374 |
| `rate-categories` `change-counts` | slow, medium, fast | 4.2e-25, 0.190, 2.137 | 5.552 |
| `rate-categories` `median-split` | slow, medium, fast | 4.4e-25, 0.513, 2.974 | 5.455 |
| `rate-categories` `k-means` | slow, medium, fast | 0.195, 2.356, 4.783 | 5.477 |

Three findings:

- **Codon position carries the rate signal; genes do not.** The codon positions span a factor of 13, the genes a factor of 1.3. Crossing the two (structural genes C, M and E against the non-structural NS1 to NS5, each split by codon position, as seven `--partition` entries) leaves the two gene groups within 10 % of each other at every codon position, so the gene split adds nothing on this genome. This is dengue, whose coding region is one polyprotein; a genome with a gene under episodic selection may differ, and `data/sc2/2844` is the test for that
- **Unsharing is right for codon positions and wrong for rate buckets.** With `--shared-params none` the codon positions fit exchangeabilities that differ by mechanism, and the differences follow the genetic code. The ratio of the C-T to the A-G exchange rate is 7.4 at first positions, where C-T is synonymous in leucine (CUA/CUG against UUA/UUG, 2 of 13 sense codons) and A-G never is (0 of 16); 0.72 at second positions, where neither is synonymous (0 of 16 and 0 of 14) and the value reflects amino-acid constraint, with A-G including the charge-preserving Lys-Arg and His-Arg swaps; and 3.42 at third positions, where both are nearly always synonymous (16 of 16 and 13 of 14) and the value is the closest estimate of the mutational spectrum. Forcing one shared W collapses all three to 3.12. The rate buckets with `--shared-params none` instead fit frequencies that are artifacts: the `fast` bucket of the `k-means` rule gives pi = [0.172, 0.395, 0.133, 0.301] against a genome composition of [0.315, 0.204, 0.261, 0.220], because conditioning on "this site changed at least 19 times" selects for states. This is evidence for the implemented defaults of `fn default_shared_params()`, which share nothing for `codon-positions` and share W and pi for `genes` and `rate-categories`
- **The rate buckets are circular, and two of the three rules are degenerate.** `change-counts` and `median-split` both define `slow` as exactly the sites with zero Fitch changes (4508 sites here, matching the 4508 constant sites that IQ-TREE reports for this alignment). A bucket of constant sites has no substitutions to explain, so its fitted rate runs to the floor: 4e-25, zero for practical purposes, and 42 % of the alignment then carries no substitution signal. Only `k-means` avoids this, because it takes its thresholds from the count distribution instead of hardcoding zero as a boundary, and pools the constant sites with enough changed sites (7662 of 10723, counts 0 to 6) to keep the rate finite at 0.195. The rates are fitted from the same change counts that defined the buckets, so the zero bucket is guaranteed by construction rather than inferred

**Recommended default for an annotated coding genome: `--partition-by codon-positions` with the default sharing (none).** The partition comes from information outside the alignment being fitted, so no quantity is used twice, unlike `rate-categories`, where the Fitch counts both define the buckets and set their rates and, when unshared, their composition. It captures the largest effect on this data, costs 24 free parameters over four partitions, and converged in 7 iterations and about 20 seconds on 1954 taxa. `rate-categories` remains the option where no annotation exists, and `k-means` is the rule to prefer there.

Two limits of the recommendation: a nucleotide GTR per codon position absorbs selection into the exchangeabilities, so a fitted W is a mixture of mutation and selection and not a mutation model (the proper fix is a codon model, [kb/issues/M-ancestral-codon-substitution-models-unimplemented.md](../issues/M-ancestral-codon-substitution-models-unimplemented.md)); and none of the comparisons above is a fit comparison, because no command reports a log-likelihood ([kb/issues/N-optimize-reports-no-log-likelihood-or-trace-log.md](../issues/N-optimize-reports-no-log-likelihood-or-trace-log.md)). The tree lengths in the table rank nothing.

### 3.4 Not implemented

- **Combined schemes.** `--partition-by` takes one scheme, so gene by codon position, or gene-group by codon position, needs hand-written `--partition` entries. The pieces exist: `owned_sites()` gives each site its CDS and `codon_position()` its codon position, and `site_ranges_from_columns()` turns any column set into entries. A combined scheme also needs a CDS-to-group mapping, which no input carries
- **Sharing by group.** `--shared-params` applies to every `infer` partition of the run. The separate W, pi and mu groups of `kb/reports/auto-partitioning.md` S7.2 are not implemented, so "share W everywhere, share pi within each codon position" cannot be expressed. `fn infer_gtr_group()` already fits one group from summed counts
- **Scheme selection.** No information criterion chooses among the schemes; see S3.3 and `kb/reports/auto-partitioning.md` S6

**Sources:** `kb/decisions/multi-alignment-partitions.md`, `packages/app-commands/src/commands/shared/partition_by.rs`, `kb/reports/auto-partitioning.md` S7.1

## 4. Model configuration (how models relate to partitions)

### 4.1 GTR family grouping specification

When multiple partitions exist, GTR parameters can be shared at three nested tiers:

- W groups (exchangeability): pool substitution counts $n_{ij}$ and time-in-state $\tau_j$ across partitions, update shared exchangeability matrix. Safe to share broadly (reflects chemistry, not selection).
- $\pi$ groups (equilibrium frequencies): pool row-summed counts and root state frequencies within groups. Share within functionally similar partitions (e.g., all 3rd codon positions).
- $\mu$ groups (overall rate): pool total counts, update shared rate. Typically per-partition (rates differ by function).

Nesting: W groups contain $\pi$ groups contain $\mu$ groups. A grouping specification defines which partitions share which parameters.

Default: share W broadly, $\pi$ within functionally similar groups, $\mu$ per-partition. Pseudocode for the generic inference function is in the auto-partitioning report.

**Sources:** `kb/reports/auto-partitioning.md` S7.2

### 4.2 Site-specific GTR as a single dense partition

Site-specific models (per-site $\pi^a$ and $\mu^a$, shared W) are one dense partition with a `SiteSpecific` model variant, not thousands of separate partitions. This preserves W sharing and traversal efficiency.

Proposed type: `enum GtrModel { Scalar(GTR), SiteSpecific(GTRSiteSpecific) }`.

v1 has no site-specific GTR code.

Design constraints: site-specific GTR is dense-only. Reject `--sparse` + `--site-specific-gtr`.

**Sources:** `kb/issues/M-gtr-site-specific-model-not-implemented.md`, `kb/reports/auto-partitioning.md` S7.3, R1, R2

### 4.3 Two-stage estimation (serialize and reload models)

Following IQ-TREE's PMSF workflow: estimate per-site preferences once on a guide tree, serialize to JSON (full model) or `.sitefreq` (IQ-TREE compatible), reload as fixed profiles for subsequent analyses. Decouples expensive estimation from repeated operations.

Two serialization proposals exist:

- **Scalar GTR:** JSON with `format: "treetime-gtr"`, fields for W (symmetric matrix), $\pi$ (frequencies), $\mu$ (rate scalar). Secondary exports: PAML triangular format, RAxML-NG model string. Enables `--model-file` and `--export-model` CLI flags.
- **Site-specific GTR:** JSON with `format: "treetime-site-specific-gtr"`, fields for W (shared matrix), $\pi$ (per-site: `n_states x seq_len`), $\mu$ (per-site: `seq_len`). Secondary export: IQ-TREE `.sitefreq`. Approximately 2 MB JSON for 30K-site nucleotide alignment.

**Sources:** `kb/proposals/model-serialization-scalar.md`, `kb/proposals/model-serialization-site-specific.md`, `kb/reports/auto-partitioning.md` S7.5

### 4.4 Edge-linked vs edge-unlinked branch lengths

Partitions sharing a tree can share branch lengths in two modes:

- Edge-proportional (scaled, RAxML-NG default): shared topology and branch lengths, per-partition rate scaler $\mu_p$. Branch length for partition $p$ on edge $e$: $b_{pe} = \mu_p \cdot b_e$.
- Edge-unlinked (independent): each partition has fully independent branch lengths. More parameters, less borrowing of signal.

v1 uses the edge-proportional mode: all partitions of optimize and timetree share the branch lengths, with one relative rate $\mu_p$ per partition (`--partition-rates`, [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)). Edge-unlinked partitions are not implemented.

**Sources:** `kb/proposals/config-file-multi-partition.md`, `kb/reports/auto-partitioning.md` S2.3, R9

## 5. Existing infrastructure

### 5.1 Partition type system

| Type | Storage | Use case |
| --- | --- | --- |
| `PartitionFitch` | Per-node state sets | Parsimony reconstruction, compression |
| `PartitionMarginalDense` | `Array2<f64>` (positions x states) per node and edge | Full probability matrices |
| `PartitionMarginalSparse` | Variable positions only | Sparse marginal reconstruction |
| `PartitionMarginalDiscrete` | Discrete trait states | Mugration |

`trait MarginalPasses` ([packages/treetime/src/partition/marginal/shared/update.rs#L10](../../packages/treetime/src/partition/marginal/shared/update.rs#L10)) holds the marginal passes that the dense, sparse and discrete partitions share. `DenseReconstruction` and `SparseReconstruction` ([packages/treetime/src/partition/marginal/reconstruction.rs](../../packages/treetime/src/partition/marginal/reconstruction.rs)) hold a partition with its GTR and its layout (alignment number and columns); `enum PartitionSet` holds the reconstructions of a run as one vector of either kind, without trait objects or locks.

**Sources:** `kb/decisions/partition-system-architecture.md`

### 5.2 Existing branches

| Branch                                     | Content                                       | Status                 |
| ------------------------------------------ | --------------------------------------------- | ---------------------- |
| `feat/abstract-partitions-enum`            | Partition abstraction via enum                | Historical, superseded |
| `feat/abstract-partitions-trait-objs`      | Partition abstraction via trait objects       | Historical, superseded |

## 6. Relationship to issue #2 (partition init duplication)

The duplication in partition init code across four commands is a symptom. The four commands do not create duplicated identical partitions. They each create independent partitions with different models, different dense/sparse decisions, and different post-init handling. The real problem is that partition creation is entangled with command orchestration instead of being a configured, composable layer.

The fix is not "extract shared init helper" but "implement partition configuration": commands receive configured partitions rather than constructing them inline. The partition entries implement this for sequence partitions:

1. `struct PartitionPlan` describes the partition to create (sites of one alignment and a model description)
2. `fn build_partition_set()` and the ancestral reconstruction plan build the partitions from the plans
3. One resolver turns the configuration sources (CLI flags, config file, app settings) into plans
4. The core pipelines receive the alignments and the plans (`struct SequenceInput`)

## 7. Cross-references

### Issues

- [N-io-multi-segment-genome-input](../issues/N-io-multi-segment-genome-input.md)
- [N-representation-infer-dense-stub](../issues/N-representation-infer-dense-stub.md)
- [kb/issues/M-gtr-site-specific-model-not-implemented.md](../issues/M-gtr-site-specific-model-not-implemented.md)
- [N-optimize-reports-no-log-likelihood-or-trace-log](../issues/N-optimize-reports-no-log-likelihood-or-trace-log.md) - blocks every fit comparison between schemes in S3.3
- [H-ancestral-on-a-time-tree-reads-years-as-substitutions](../issues/H-ancestral-on-a-time-tree-reads-years-as-substitutions.md) - a fitted model cannot be carried to a downstream command with its rate
- [M-optimize-loop-stops-on-decreased-likelihood-with-long-branches](../issues/M-optimize-loop-stops-on-decreased-likelihood-with-long-branches.md) - triggered by `rate-categories` in the measurements of S3.3

### Proposals

- [config-file-multi-partition](../proposals/config-file-multi-partition.md)
- [model-serialization-scalar](../proposals/model-serialization-scalar.md)
- [model-serialization-site-specific](../proposals/model-serialization-site-specific.md)

### Decisions

- [partition-system-architecture](../decisions/partition-system-architecture.md)
- [multi-alignment-partitions](../decisions/multi-alignment-partitions.md)

### Reports

- [auto-partitioning](auto-partitioning.md)

### Features

- [representation](../features/representation.md) - automatic codon-position partitioning `[ ]`
