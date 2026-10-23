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

### 2.4 `--segment` flag (multi-segment genomes)

For segmented viruses (influenza: 8 segments, same alphabet, independent evolution), a `--segment` flag would load segment-aware FASTA and create one partition per segment. Segments that share one tree run as one `--alignment` file per segment, or as `--partition` entries of one file (S2.3); segments with different trees remain open.

v0 handles segments via concatenation.

**Sources:** `kb/issues/N-io-multi-segment-genome-input.md`

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

## 3. Automatic partitioning (heuristic derivation)

Beyond user specification, partitions can be derived automatically from alignment properties. These methods take one alignment and produce a `position -> partition_index` map.

### 3.1 Rate-based partitioning (mutation count quantiles)

Run Fitch parsimony on the full alignment. Count mutations per site. Sort sites by mutation count. Split into $q$ quantile bins. Sites in each bin form one partition. Related to k-means site-rate clustering (Frandsen et al. 2015).

This creates partitions where within-partition rate variation is minimized, improving GTR model fit per partition.

**Sources:** `kb/reports/auto-partitioning.md` S7.1

### 3.2 Codon-position partitioning (from GFF annotation)

Parse GFF annotations to identify CDS regions. Split protein-coding positions by codon position (1st, 2nd, 3rd). Non-coding positions form a 4th partition. Can further isolate specific genes (e.g., Spike protein) into their own 1st/2nd/3rd codon-position partitions.

Biologically motivated: 3rd codon positions evolve faster (wobble), 1st/2nd positions are more constrained. Mixing them in one GTR model distorts both rate and frequency estimates.

Feature inventory status: `[ ]` not implemented.

**Sources:** `kb/reports/auto-partitioning.md` S7.1, `kb/features/representation.md`

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
